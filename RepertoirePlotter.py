import matplotlib.pyplot as plt
import numpy as np
import pyphi
import os

class RepertoirePlotter:
    '''This is a class providing static methods to automate plotting of
    PyPhi's repertoire probability distributions.
    '''

    @staticmethod
    def plot(repertoire, xlabs, title, path, fname):
        '''Common function shared by both plot_from_tpm() and plot_core().

        Parameters:
            repertoire (np.ndarray): 1D vector with probabilities.
            xlabs (list[str]): Labels used for the x-axis.
            title (str): Title used in the plot.
            path (str): Directory where plot will be saved.
            fname (str): Filename and figure title.

        Returns:
            None
        '''

        fig, ax = plt.subplots()
        rects = ax.bar(range(len(repertoire)), repertoire)
        yticks = np.arange(0, 1, .25)
        plt.xticks(range(len(repertoire)))
        plt.yticks(yticks)
        ax.set_xticklabels(xlabs, size='xx-large', rotation=0)
        ax.set_yticklabels(yticks, size='xx-large')
        ax.set_title(title)

        # annotate each bar with its probability value
        for rect in rects:
            height = rect.get_height()
            ax.annotate('{}'.format(np.round(height, 4)),
                        xy=(rect.get_x() + rect.get_width() / 2, height),
                        xytext=(0, 3),  # 3 points vertical offset
                        textcoords="offset points",
                        ha='center',
                        va='bottom',
                        size='xx-large')

        # display plot and also save to file
        os.makedirs(path, exist_ok = True)
        plt.savefig(path + "/" + fname + ".svg", transparent=False)
        plt.show()

    @staticmethod
    def plot_from_tpm(tpm, given_state, other_time_states, xlabs, path, fname, effect_purview=True):
        '''Function to plot candidate cause and effect repertoires given a
        system's TPM, the selected mechanism, a state and candidate
        purview.

        Parameters:
            tpm (np.ndarray): The transition probability matrix of the
                              network, in state-by-state form.
            given_state (list[int]): row or column indices to slice the TPM.
            other_time_states (list[list[int]]): how to group states at
                                                 time A in Pr(S_A|S_B).
            xlabs (list[str]): labels used for the x-axis.
            path (str): directory where plot will be saved.
            fname (str): filename and figure title.
            effect_purview (bool): whether purview is effect or cause.

        Returns:
            None

        Example:
            # effect of [289] over [288, 289]
            RepertoirePlotter.plot_from_tpm(TPM_max_phi_coarse_grained,
                                            given_state=[0,1],
                                            other_time_states=[[0],[1],[2],[3]],
                                            xlabs=["OFF,OFF", "ON,OFF", "OFF,ON", "ON,ON"],
                                            path="figures/fig4/",
                                            fname=r'$P(S_{system_{t+1}} | S_{289_t} = OFF)$',
                                            effect_purview=True)
        '''

        # slice conditional distribution depending on current mechanism state
        if effect_purview:
            repertoire = np.sum(tpm[given_state, :], axis=0)
        else:
            repertoire = np.sum(tpm[:, given_state], axis=1)

        repertoire = repertoire / np.sum(repertoire)

        # aggregate "other-time" elements according to purview
        repertoire2 = []
        for micros in other_time_states:
            repertoire2.append(np.sum(repertoire[micros]))

        repertoire = repertoire2
        RepertoirePlotter.plot(repertoire,
                               xlabs,
                               fname + " = " + str(np.round(repertoire, 4)),
                               path,
                               fname)

    @staticmethod
    def xlabs_from_repertoire(repertoire):
        '''Ugly helper function to semi-automate x-labels generation.'''

        if repertoire.size == 2:
            return ["OFF", "ON"]
        elif repertoire.size == 4:
            return ["OFF,OFF", "ON,OFF", "OFF,ON", "ON,ON"]
        else:
            raise Exception('not implemented')

    @staticmethod
    def plot_core(concept, path, MIP=False):
        '''Function to plot core cause and effect repertoires after running
        pyphi.compute.sia().

        Parameters:
            concept (pyphi.models.mechanism.Concept): concept/mechanism to plot.
            path (str): directory where plot will be saved.
            MIP (bool): whether the repertoire for the Minimum Information
                        Partition should be plotted instead.

        Returns:
            None

        Example:
            network = pyphi.Network(TPM_max_phi_coarse_grained,
                        num_states_per_node=num_states_l[0])
            state = (0,0)
            subsystem = pyphi.Subsystem(network, state)
            sia = pyphi.compute.sia(subsystem)

            # Just one concept/mechanism (sia.ces[0])
            concept=sia.ces[0]
            # "normal repertoire"
            RepertoirePlotter.plot_core(concept, "save/at/path/", MIP=False)
            # "partitioned repertoire"
            RepertoirePlotter.plot_core(concept, "save/at/path/", MIP=True)
        '''
        if not MIP:
            repertoire = concept.cause.repertoire
        else:
            repertoire = concept.cause.partitioned_repertoire
        xlabs=RepertoirePlotter.xlabs_from_repertoire(repertoire)
        RepertoirePlotter.plot(repertoire.reshape([repertoire.size], order='F'),
                               xlabs,
                               str(concept.cause.phi) + "\n" + str(concept.cause.mip),
                               path,
                               fname='cause' if not MIP else 'cause (MIP)')
        if not MIP:
            repertoire = concept.effect.repertoire
        else:
            repertoire = concept.effect.partitioned_repertoire
        xlabs=RepertoirePlotter.xlabs_from_repertoire(repertoire)
        RepertoirePlotter.plot(repertoire.reshape([repertoire.size], order='F'),
                                          xlabs,
                                          str(concept.cause.phi) + "\n" + str(concept.cause.mip),
                                          path,
                                          fname='effect' if not MIP else 'effect (MIP)')
