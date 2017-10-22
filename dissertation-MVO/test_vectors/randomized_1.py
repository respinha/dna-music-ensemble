from ensemble import Composer, FileWriter
from algorithms import PitchAlgorithm, DurationsAlgorithm, DynamicsAlgorithm, ClusteringAlgorithm
from music21 import scale
import numpy as np


def run():

    alphabet = ['a','c','g','t','-']

    np.random.choice(alphabet, 50, p=[0.5, 0.1, 0.1, 0.3])
    msa = np.zeros((2,200), dtype="S1")

    msa[0] = np.zeros(200, dtype="S1")

    msa[0][0:50] = 'a'
    msa[0][50:75] = 'c'
    msa[0][75:150] = '-'
    msa[0][150:175] = 'g'
    msa[0][150:200] = '-'

    print('MSA 0', msa[0])
    msa[1][0:25] = 'a'
    msa[1][25:100] = 'c'
    msa[1][100:175] = 'g'
    msa[1][175:200] = 'a'

    print('MSA 1', msa[1])

    pitch_algo = PitchAlgorithm(PitchAlgorithm.WORD_DISTANCES, scale=scale.MajorScale().getPitches(), n_nucleotides=1)
    durations_algo = DurationsAlgorithm(DurationsAlgorithm.FREQUENCIES_DYNAMIC, window_size=10, window_duration=10, n_nucleotides=1)
    dynamics_algo = DynamicsAlgorithm(DynamicsAlgorithm.SHANNON_INDEX, window_size=10, gap_window_threshold=0.5, gap_column_threshold=0.7, criteria='local', levels=7)
    instruments_algo = ClusteringAlgorithm('kmeans')


    composer = Composer(instruments_algo, pitch_algo, durations_algo, dynamics_algo, input_type='array',
                        alignment=msa)

    sequence_names = np.array([str(i) for i in range(0, len(msa))])

    i = 0

    scores = composer.gen_numerical_vectors(k=2, piece_length=20)
    for score in scores:
        fw = FileWriter(score, sequence_names)
        fname = 'demo_' + str(i)

        print('fname', fname)
        fw.write(midi=fname, audio=fname, display=False)
        i += 1

if __name__ == '__main__':
    run()