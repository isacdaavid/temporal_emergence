import matplotlib.pyplot as plt
import numpy as np
import pyphi
import os

class RepertoirePlotter:

    @staticmethod
    def plot_from_tpm(tpm, given_state, other_time_states, xlabs, path, fname, effect_purview=True):
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

        # plot
        fig, ax = plt.subplots()
        ax.bar(range(len(repertoire)), repertoire)
        yticks = np.arange(0, 1, .25)
        plt.xticks(range(len(repertoire)))
        plt.yticks(yticks)
        ax.set_xticklabels(xlabs, size='xx-large', rotation=0)
        ax.set_yticklabels(yticks, size='xx-large')
        ax.set_title(fname + " = " + str(np.round(repertoire, 4)))
        os.makedirs(path, exist_ok = True)
        plt.savefig(path + fname + ".svg", transparent=False)
        plt.show()

    @staticmethod
    def plot(repertoire, xlabs, phi, mip, path, fname):
        fig, ax = plt.subplots()
        rects = ax.bar(range(len(repertoire)), repertoire)
        yticks = np.arange(0, 1, .25)
        plt.xticks(range(len(repertoire)))
        plt.yticks(yticks)
        ax.set_xticklabels(xlabs, size='xx-large', rotation=0)
        ax.set_yticklabels(yticks, size='xx-large')
        ax.set_title(f"φ = {str(phi)}\n{str(mip)}")

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
    def xlabs_from_repertoire(repertoire):
        if repertoire.size == 2:
            return ["OFF", "ON"]
        elif repertoire.size == 4:
            return ["OFF,OFF", "ON,OFF", "OFF,ON", "ON,ON"]
        else:
            raise Exception('not implemented')

    @staticmethod
    def plot_core(concept, path, MIP=False):
        if not MIP:
            repertoire = concept.cause.repertoire
        else:
            repertoire = concept.cause.partitioned_repertoire
        xlabs=RepertoirePlotter.xlabs_from_repertoire(repertoire)
        RepertoirePlotter.plot(repertoire.reshape([repertoire.size], order='F'),
                               xlabs,
                               concept.cause.phi,
                               concept.cause.mip,
                               path,
                               fname='cause' if MIP == False else 'cause (MIP)')
        if not MIP:
            repertoire = concept.effect.repertoire
        else:
            repertoire = concept.effect.partitioned_repertoire
        xlabs=RepertoirePlotter.xlabs_from_repertoire(repertoire)
        RepertoirePlotter.plot(repertoire.reshape([repertoire.size], order='F'),
                                          xlabs,
                                          concept.effect.phi,
                                          concept.effect.mip,
                                          path,
                                          fname='effect' if MIP == False else 'effect (MIP)')
