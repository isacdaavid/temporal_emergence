import numpy as np 
import matplotlib.pyplot as plt
from numpy.lib.arraysetops import isin
import pandas as pd
import random
import math
import pyphi # needs nonbinary install
pyphi.config.PARTITION_TYPE = 'ALL'
# pyphi.config.MEASURE = 'AID'
# pyphi.config.USE_SMALL_PHI_DIFFERENCE_FOR_CES_DISTANCE = True
pyphi.config.ASSUME_CUTS_CANNOT_CREATE_NEW_CONCEPTS = True
pyphi.config.WELCOME_OFF = True

class Neuron:

    @staticmethod
    def binarise_spiketrain(a,S):
        """
        Binarises a spike-train with bins of size S. 
        - a is a list of times when the neuron fired. 
        """
        a_states = np.zeros(math.ceil(max(a)/S))
        for spike_t in a:
            index = int(spike_t / S)
            a_states[index] = 1
        return a_states

class TPMMaker:

    @staticmethod    
    def get_TPM_index(state):  # rename this
        """
        Given the state of a set of neurons, get the
        index into the TPM that the state corresponds to.

        Example: 
            state = [[1,0,0]
                    [0,1,1]]
            Neuron A is in state 100 and B is in state 011
            State A corresponds to the decimal 4
            State B corresponds to the decimal 3
            The index of A, the first listed state, varies 'fastest' in the TPM according to the PyPhi convention, 
            https://pyphi.readthedocs.io/en/latest/conventions.html.
            So, the indices, as tuples, will be (A0,B0),(A1,B0),(A2,B0),(A3,B0),(A4,B0),...,(A0,B1),(A1,B1)... 
            To convert these tuples to ints, we need A + B*8 = 4 + 3*8 = 28

            In general, where each node has n possible states, this becomes:
                - A*n^0 + B*n^1 + C*n^2 + ... where A, B, C are the states of each node,
                in order from first to last in the state input. 
                - Note that for the binary arrays accepted by this function, 
                n = 2^(m), where m is the size of the binary array that describes the state of a node. 
        """
        n = 2**(state.shape[1])
        index = 0
        for i in range(state.shape[0]):
            dec_node = int("".join(str(int(j)) for j in state[i,:]), 2)
            index += dec_node * n ** i
        return index
    @staticmethod
    def get_num_state_occurrences(spiketrains, S, K, skipby):
        """Gets the number of occurrences of each state in the TPM, if we sample
        all occurrences of each state when creating the TPM (which we don't)
        """
        binaryneurons = TPMMaker.get_binarised_trains(spiketrains, S)
        assert K >= 1
        assert binaryneurons.shape[0] >= 1

        size = (2**K)**binaryneurons.shape[0]
        
        # initialise num_transition 
        num_transitions = np.zeros(size)

        # get a randomly ordered list of indices at which to look at transitions
        # start at K-1 because our state at time i looks BACK to i-1, i-2,.. to build the rest of state
        rand_indices = np.array(list(range(K-1, binaryneurons.shape[1] - skipby)))
        np.random.shuffle(rand_indices)

        for i in rand_indices:
            curr_state = binaryneurons[:,i-(K-1):(i+1)]
            i_c = TPMMaker.get_TPM_index(curr_state)
            num_transitions[i_c] += 1
        
        return num_transitions
    @staticmethod
    def get_TPM_nonbinary(binaryneurons, K, skipby, required_obs):
        """Given an array of binarised neuron spike-trains
        and a K value for how many time-steps to include in a single state, 
        get the TPM of the system. 
            - Skipby controls where the future state starts: given that the current state
            starts at T, future state starts at T+skipby
        Example: 
            if K = 3, then find the TPM that describes 
            the transition probability of System[t-2,t-1,t] --> System[t+1,t+2,t+3]

        Returns:
            - A TPM of the system in state-state mode (TODO: conventions?)
            - A matrix A with the same dimensions of TPM, where A[i,j] is the 
            number of transitions that were used to calculate the value of TPM[i,j].
        
        """
        assert K >= 1
        assert binaryneurons.shape[0] >= 1

        size = (2**K)**binaryneurons.shape[0]
        
        # initialise TPM and num_transitions arrays
        TPM = np.zeros((size, size))
        num_transitions = np.zeros((size, size))

        # get a randomly ordered list of indices at which to look at transitions
        # start at K-1 because our state at time i looks BACK to i-1, i-2,.. to build the rest of state
        rand_indices = np.array(list(range(K-1, binaryneurons.shape[1] - skipby)))
        np.random.shuffle(rand_indices)
        for i in rand_indices:
            curr_state = binaryneurons[:,i-(K-1):(i+1)]
            i_c = TPMMaker.get_TPM_index(curr_state)
            total = sum(num_transitions[i_c,:])
            if total >= required_obs:   # don't add this observation if we already have enough
                continue

            future_state = binaryneurons[:,i-(K-1) + skipby:(i+1) + skipby]   # Ugly indexing 
            i_f = TPMMaker.get_TPM_index(future_state)

            num_transitions[i_c, i_f] += 1
        
        for j in range(num_transitions.shape[0]):
            total = sum(num_transitions[j,:])
            if total < required_obs:
                raise ValueError("State with index " + str(j) + \
                " was observed in the data fewer than " + str(required_obs) + " times, (" + str(total) + " times only).")
            
            TPM[j,:] = num_transitions[j,:] / total
        
        return TPM, num_transitions
        

    @staticmethod
    def get_TPM_nonbinary_nonrandom(binaryneurons, K, skipby, required_obs, startindex):
        """Given an array of binarised neuron spike-trains
        and a K value for how many time-steps to include in a single state, 
        get the TPM of the system. 
            - Skipby controls where the future state starts: given that the current state
            starts at T, future state starts at T+skipby
        Example: 
            if K = 3, then find the TPM that describes 
            the transition probability of System[t-2,t-1,t] --> System[t+1,t+2,t+3]

        Returns:
            - A TPM of the system in state-state mode (TODO: conventions?)
            - A matrix A with the same dimensions of TPM, where A[i,j] is the 
            number of transitions that were used to calculate the value of TPM[i,j].
        
        - always samples in order, starting from startindex
        """
        assert K >= 1
        assert binaryneurons.shape[0] >= 1

        size = (2**K)**binaryneurons.shape[0]
        
        # initialise TPM and num_transitions arrays
        TPM = np.zeros((size, size))
        num_transitions = np.zeros((size, size))

        # get an ordered list of indices at which to look at transitions
        # start at K-1 because our state at time i looks BACK to i-1, i-2,.. to build the rest of state
        # but here we actually want to start at K-1 + startindex to shift to startindex 
        rand_indices = np.array(list(range(K-1+startindex, binaryneurons.shape[1] - skipby, 2)))
        #print(rand_indices)
        for i in rand_indices:
            curr_state = binaryneurons[:,i-(K-1):(i+1)]
            i_c = TPMMaker.get_TPM_index(curr_state)
            total = sum(num_transitions[i_c,:])
            if total >= required_obs:   # don't add this observation if we already have enough
                continue

            future_state = binaryneurons[:,i-(K-1) + skipby:(i+1) + skipby]   # Ugly indexing 
            i_f = TPMMaker.get_TPM_index(future_state)

            num_transitions[i_c, i_f] += 1
        
        for j in range(num_transitions.shape[0]):
            total = sum(num_transitions[j,:])
            if total < required_obs:
                raise ValueError("State with index " + str(j) + \
                " was observed in the data fewer than " + str(required_obs) + " times, (" + str(total) + " times only).")
            
            TPM[j,:] = num_transitions[j,:] / total
        
        return TPM, num_transitions

    @staticmethod
    def get_binarised_trains(spiketrains,S):
        # get the binarised spike trains for each neuron from a train of float spikes
        # create a multidimensional array of the spiketrains by not considering further than the shortest train
        binarised_trains = [[] for _ in range(len(spiketrains))]
        min_length = math.inf
        for i in range(len(spiketrains)):
            binarised_trains[i] = Neuron.binarise_spiketrain(spiketrains[i], S)
            min_length = min(min_length, len(binarised_trains[i]))
        
        return np.array([binarised_trains[i][0:min_length] for i in range(len(binarised_trains))])   # probably slow? 


    @staticmethod
    def TPM_from_spiketrains(spiketrains, S, K, skip, required_obs):
        binarised_trains = TPMMaker.get_binarised_trains(spiketrains, S)
        return TPMMaker.get_TPM_nonbinary(binarised_trains, K, skip, required_obs)


class CoarseGrainer:

    @staticmethod
    def coarse_grain_nonbinary_TPM(TPM, state_map, num_states_per_elem):
        """
        TODO: add checks to inputs

        Turns a nonbinary TPM into a binary TPM by coarse-graining each nonbinary element
        according to a grouping.
            - grouping only groups states of individual elements.
            - grouping example for a TPM of 2 elements each w 4 states:
            [[[0, 1], [2, 3]], [[0], [1, 2], [3]]]
                - The first element's 0 and 1 states are grouped into a state, and 2,3 into another.
                - The second element's 0 state stays the same, 1,2 are grouped and 3 stays the same.

        state_map:
            - a dictionary that takes in the index of a macro state, obtained from num_states_per_elem
            and returns a list of micro states to the TPM that make up the macro state. 

        For example, for a TPM of 2 elements, each with 4 states,
        the input TPM will have 4*4 = 16 states. We can coarse grain
        such that only the first three states of each element will be grouped into an OFF state,
        and the last element will be grouped into ON. 

        However, we don't have to coarsegrain to binary. We might say that the first state
        of each element will map to OFF, the second and third will map to FIRING, 
        fourth to BURSTING. 
        """
        num_states = 1
        for i in num_states_per_elem:
            num_states *= i

        new_TPM = np.zeros((num_states, num_states))
        for i in range(num_states):
            input_micro_indices = state_map[i]
            for j in range(num_states):
                output_micro_indices = state_map[j]
                
                total_prob = 0
                """
                for output_state in output_micro_indices:
                    avg_for_output_state = 0
                    for input_state in input_micro_indices:
                        avg_for_output_state += TPM[input_state, output_state]

                    avg_for_output_state = avg_for_output_state / len(input_micro_indices)
                    total_prob += avg_for_output_state
                """
                for input_state in input_micro_indices:
                    total_for_input_state = 0
                    for output_state in output_micro_indices:
                        total_for_input_state += TPM[input_state, output_state]
                    
                    total_prob += total_for_input_state
                new_TPM[i,j] = total_prob / len(input_micro_indices)
                
        return new_TPM

    @staticmethod
    def get_state_map(coarse_grain):
        """Get the map from each micro state to its macro state
        IMPORTANT: Works only for coarse grainings of pairs of neurons! 
        """
        assert len(coarse_grain) == 2, "Require exactly 2 elements"
        
        num_micro_states_per_elem = []
        num_macro_states = 1
        for elem in coarse_grain:
            num_macro_states *= len(elem)
            num_micro_states_per_elem.append(sum(len(x) for x in elem))

        
        new_coarse_grain = []
        for elem in coarse_grain:
            states = []
            for s in range(len(elem)):
                states.append((elem[s], s))
            new_coarse_grain.append(states)

        # for every element pair
        state_map = {}
        for q in range(len(new_coarse_grain[1])):
            for p in range(len(new_coarse_grain[0])):
                states_p = list(new_coarse_grain[0][p][0])
                states_q = [num_micro_states_per_elem[0] * e for e in new_coarse_grain[1][q][0]]
                
                micro_system_states = np.add.outer(states_p, states_q).ravel()
                macro_state = p + q * len(new_coarse_grain[0])
                state_map[macro_state] = micro_system_states.tolist()
        return state_map

    @staticmethod
    def get_state_maps(element_coarse_grainings):
        """For a 2 element system, given coarse graining options for an element,
           get the state maps for each coarse graining combination"""
        states = []
        num_states_l = []
        for e1 in element_coarse_grainings:
            for e2 in element_coarse_grainings:
                c_state_map = CoarseGrainer.get_state_map([e1, e2])
                c_num_states = [len(e) for e in [e1,e2]]

                states.append(c_state_map)
                num_states_l.append(c_num_states)
        return states, num_states_l

class PhiCalculator:

    @staticmethod
    def get_micro_phis(TPM, verbose=True):
        """
        Gets the state phis for a TPM with 2 elements, each 4 states.
        """
        network = pyphi.Network(
        (TPM),
        num_states_per_node=[4,4]
        )
        phis = []
        states = Helpers.get_nary_states(2,4)   # not the TPM order! 
        for state in states:
            subsystem = pyphi.Subsystem(network, state)
            sia = pyphi.compute.sia(subsystem)
            if state == (0,1):
                if verbose:
                    print(sia.ces)
                    print(sia.partitioned_ces)
                    print(sia.cut)
            phis.append(sia.phi)
        if verbose:
            print(phis)
        
        return phis
    
    @staticmethod
    def get_micro_average_phi(TPM,verbose=True):
        phis = PhiCalculator.get_micro_phis(TPM, verbose)
        return sum(phis) / len(phis)

    @staticmethod
    def get_micro_weighted_average_phi(TPM,occurrences,verbose=True):
        """
        Weighted average of phis, weighting state phis by the occurrence probabilities
        of each state. 
        """
        phis = PhiCalculator.get_micro_phis(TPM, verbose)
        print(phis)
        total_occs = sum(occurrences)
        weights = [o / total_occs for o in occurrences]
        return sum([phis[i] * weights[i] for i in range(len(phis))]), weights

    @staticmethod
    def get_macro_phis(micro_TPM, verbose=True, state_map=None, num_states_per_elem=None):

        if state_map == None or num_states_per_elem == None:   # have a default state map
            state_map = {0: [0,1,2, 4,5,6, 8,9,10], 1: [3,7,11], 2: [12,13,14], 3:[15]} 
            num_states_per_elem = [2,2]
        
        macro_TPM = CoarseGrainer.coarse_grain_nonbinary_TPM(micro_TPM, state_map, num_states_per_elem)
        network = pyphi.Network(
        macro_TPM,
        num_states_per_node=num_states_per_elem
        )

        states = Helpers.get_system_states(num_states_per_elem)
        phis = []
        for state in states:
            subsystem = pyphi.Subsystem(network, state)
            sia = pyphi.compute.sia(subsystem)
            phis.append(sia.phi)
        
        if verbose:
            print(sia.ces)
            print(sia.partitioned_ces)
            print(sia.cut)
        
        return phis

    @staticmethod
    def get_macro_average_phi(micro_TPM, verbose=True, state_map=None, num_states_per_elem=None):
        phis = PhiCalculator.get_macro_phis(micro_TPM, verbose, state_map, num_states_per_elem)
        return sum(phis) / len(phis)

    @staticmethod
    def get_macro_weighted_average_phi(micro_TPM, occurrences, verbose=True, state_map=None, num_states_per_elem=None):
        phis = PhiCalculator.get_macro_phis(micro_TPM, verbose, state_map, num_states_per_elem)
        # TODO get rid of repeated code
        if state_map == None or num_states_per_elem == None:   # have a default state map
            state_map = {0: [0,1,2, 4,5,6, 8,9,10], 1: [3,7,11], 2: [12,13,14], 3:[15]} 
            num_states_per_elem = [2,2]

        macro_occurrences = [sum([occurrences[i] for i in state_map[key]]) for key in state_map]
        total_occurrences = sum(macro_occurrences)
        weights = [m / total_occurrences for m in macro_occurrences]
        weighted_phis = [phis[i]  * weights[i] for i in range(len(phis))]

        return sum(weighted_phis) / len(weighted_phis)

    @staticmethod
    def all_coarsegrains_get_macro_average_phi(micro_TPM, verbose=True):
        
        # ways to coarse grain each element
        element_coarse_grainings = [[[0], [1,2,3]], [[0,1,2],[3]], [[0], [1,2], [3]], [[0], [1], [2], [3]]]
        states, num_states_l = CoarseGrainer.get_state_maps(element_coarse_grainings)
            
        phis = []
        for i in range(len(states)):
            state_map, num_states = states[i], num_states_l[i]
            if verbose:
                print(state_map, num_states)
                print("\n")
            average_phi = PhiCalculator.get_macro_average_phi(micro_TPM, verbose=verbose, state_map=state_map, num_states_per_elem=num_states)
            phis.append(average_phi)
        
        return phis
    
    @staticmethod
    def all_coarsegrains_get_macro_weighted_average_phi(micro_TPM, occurrences, verbose=True):
        # ways to coarse grain each element
        element_coarse_grainings = [[[0], [1,2,3]], [[0,1,2],[3]], [[0], [1,2], [3]], [[0], [1], [2], [3]]]
        states, num_states_l = CoarseGrainer.get_state_maps(element_coarse_grainings)
            
        phis = []
        for i in range(len(states)):
            state_map, num_states = states[i], num_states_l[i]
            if verbose:
                print(state_map, num_states)
                print("\n")
            average_phi = PhiCalculator.get_macro_weighted_average_phi(micro_TPM, occurrences, verbose=verbose, state_map=state_map, num_states_per_elem=num_states)
            phis.append(average_phi)
        return phis


class Helpers:

    @staticmethod
    def get_bin_states(l):
        states = []
        for i in range(2**l):
            state = list(bin(i)[2:])
            while len(state) < l:
                state.insert(0,'0')
            state = tuple([int(c) for c in state])
            states.append(state)
        return states

    # https://stackoverflow.com/a/28666223
    @staticmethod
    def numberToBase(n, b):
        if n == 0:
            return [0]
        digits = []
        while n:
            digits.append(int(n % b))
            n //= b
        return digits[::-1]

    @staticmethod
    def get_nary_states(n,b):
        states = []
        for i in range(b**n):
            state = Helpers.numberToBase(i,b)
            while len(state) < n:
                state.insert(0,0)
            state = tuple(state)
            states.append(state)
        return states
    
    @staticmethod
    def get_system_states(state_nums):
        elem_states = [list(range(i)) for i in state_nums]
        states = elem_states[0]
        for i in range(1,len(elem_states)):
            states = Helpers.outer(states, elem_states[i], Helpers.join_ints)
        return states
    
    @staticmethod
    def outer(l1, l2, f):
        res = []
        for i in l1:
            for j in l2:
                res.append(f(i,j))

        return res
    
    @staticmethod
    def join_ints(i1, i2):
        if isinstance(i1, tuple):
            return (*i1,i2)
        else:
            return (i1, i2)
    
    @staticmethod
    def optional_asterisk(tup):
        """If tup is a tuple, then unpack it with asterisk, but otherwise 
        just leave it be
        """
        
        
class DataGenerator:
    def __init__(self, TPM, base):
        self.TPM = TPM
        # base is the number of possible states in each state of the TPM
        self.base = base
        self.num_nodes = int(math.log(self.TPM.shape[0], self.base))

    def generate_timeseries(self, iters):
        bits_in_state = int(math.log2(self.base))
        data = np.zeros((self.num_nodes, iters*bits_in_state))
        curr_state = 0
        for i in range(0, iters*bits_in_state, bits_in_state):
            # see https://stackoverflow.com/a/41852266
            curr_state = np.random.choice(np.arange(len(self.TPM[curr_state])), p=self.TPM[curr_state])

            # get the system state in array form, with the state of each node 
            state_array = Helpers.numberToBase(curr_state, self.base)
            while len(state_array) < self.num_nodes:
                state_array.insert(0,0)
            state_array = state_array[::-1]

            for j in range(len(state_array)):   # for each node
                # we need to set multiple values of the timeseries, because each iteration in the nonbinary TPM 
                # is actually multiple binary steps, the nonbinary representation 'bunches up' binary data.
                binary_state = Helpers.numberToBase(state_array[j], 2)
                while len(binary_state) < bits_in_state:
                    binary_state.insert(0,0)
                for k in range(bits_in_state):
                    data[j, i+k] = binary_state[k]
        return data

def get_phis(r, t, num_transitions, infolder, outfolder):
    ### LOAD DATASET ###
    print("get_phis")
    i_sec = np.loadtxt(infolder + "/cell" + str(r) + ".txt") / 1000   # divide through as they are loaded in miliseconds
    j_sec = np.loadtxt(infolder + "/cell" + str(t) + ".txt") / 1000
    cluster = np.array([i_sec, j_sec])
    print("after load dataset")
    ### COMPUTE PHIS ###
    NUM_COARSE_GRAININGS = 16
    NUM_BITS = 2
    skips = list(range(2,11,2))

    max_binsize = 0.02  # 20 ms bins
    min_binsize = 0.0029 # skip 1ms bins  -   never work and are very slow to compute
    num_binsizes = 9
    binsizes = np.linspace(min_binsize, max_binsize, num_binsizes)


    micro_phis = np.zeros((len(binsizes), len(skips)))
    macro_phis = [np.zeros((len(binsizes), len(skips),NUM_COARSE_GRAININGS))]

    for i in range(len(binsizes)):
        binsize = binsizes[i]
        for j in range(len(skips)):
            skip = skips[j]

            try:
                TPM,_ = TPMMaker.TPM_from_spiketrains(cluster,binsize,NUM_BITS,skip,num_transitions)
                tpmname = "micro_" + str(i) + "_" + str(j) + "_occs_" + str(num_transitions) + "_bin_"+str(binsize)+"_skip_"+str(skip)+".csv" 
                np.savetxt(outfolder+"/"+tpmname, TPM)
                success = True
            except:
                success = False
            
            if success:
                micro_phis[i,j] = PhiCalculator.get_micro_average_phi(TPM, verbose=False)
                all_coarse_macros = PhiCalculator.all_coarsegrains_get_macro_average_phi(TPM, verbose=False)
                macro_phis[i,j] = all_coarse_macros
                #macro_phis[i,j] = np.nanmax(all_coarse_macros)
            
            else:
                micro_phis[i,j] = None
                macro_phis[i,j] = [None for i in range(NUM_COARSE_GRAININGS)]
    
    micro_phis = np.array(micro_phis, dtype=np.float64)
    macro_phis = np.array(macro_phis, dtype=np.float64)

    np.save(outfolder + "/micro_" + str(r) + "_" + str(t), micro_phis)
    np.save(outfolder + "/macro_" + str(r) + "_" + str(t), macro_phis)
    #max_micro = np.nanmax(micro_phis)
    #max_macro = np.nanmax(macro_phis)
    #macro_win = True if max_macro > max_micro else False

    return (micro_phis, macro_phis, (r,t))



if __name__ == "__main__":
    
    ## LOAD DATA ## 
    n_143 = np.loadtxt("GLMCC/Cori_2016-12-14_probe1/cell143.txt") / 1000
    n_168 = np.loadtxt("GLMCC/Cori_2016-12-14_probe1/cell168.txt") / 1000
    a = Neuron.binarise_spiketrain(n_143, 0.002)
    b = Neuron.binarise_spiketrain(n_168, 0.002)
    print(len(a))
    print(len(b))
    
    n_168 = np.loadtxt("GLMCC/Cori_2016-12-14_probe1/cell168.txt") / 1000
    cluster_143_168 = np.array([n_168,n_143])

    NUM_BITS = 2
    skips = [2]#list(range(2,11,2))

    max_binsize = 0.02  # 50 ms bins
    min_binsize = 0.001 # 1ms bins  -   probably won't work
    num_binsizes = 10
    binsizes = [0.002]#np.linspace(min_binsize, max_binsize, num_binsizes)


    num_transitions = 200
    micro_phis = np.zeros((len(binsizes), len(skips)))
    macro_phis = np.zeros((len(binsizes), len(skips)))

    for i in range(len(binsizes)):
        binsize = binsizes[i]
        for j in range(len(skips)):
            skip = skips[j]

            try:
                TPM,_ = TPMMaker.TPM_from_spiketrains(cluster_143_168,binsize,NUM_BITS,skip,num_transitions)
                occurrences = TPMMaker.get_num_state_occurrences(cluster_143_168,binsize,NUM_BITS,skip)
                #tpmname = "micro_143_168_bin_"+str(binsize)+"_skip_"+str(skip)+".csv" 
                #np.savetxt("TPMs/"+tpmname, TPM)
                print(occurrences)

                success = True
            except ValueError:
                success = False
                print("Failed for binsize: " + str(binsize) + " and skip: " + str(skip))
            
            if success:
                #micro_phis[i,j] = PhiCalculator.get_micro_average_phi(TPM, verbose=False)
                #macro_phis[i,j] = PhiCalculator.get_macro_average_phi(TPM, verbose=False)
                micro = PhiCalculator.get_micro_weighted_average_phi(TPM, occurrences, verbose=False)
                macros_weighted = PhiCalculator.all_coarsegrains_get_macro_weighted_average_phi(TPM, occurrences, verbose=False)
                macros = PhiCalculator.all_coarsegrains_get_macro_average_phi(TPM, verbose=False)
                print(macros)
                print(macros_weighted)
                print("Success for binsize: " + str(binsize) + " and skip: " + str(skip))
            
            else:
                micro_phis[i,j] = None
                macro_phis[i,j] = None
    
    """
    print(list(reversed(Helpers.get_system_states([2,2]))))
    states = reversed(Helpers.get_nary_states(2, 2)) # not the TPM order! 
    print(list(states))
    """
    #print(CoarseGrainer.get_state_map([[[0], [1,2,3]], [[0], [1,2], [3]]]))
