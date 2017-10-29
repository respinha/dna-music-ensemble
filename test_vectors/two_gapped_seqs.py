from config import GLOBALS
from ensemble import Composer, FileWriter
from algorithms import PitchAlgorithm, DurationsAlgorithm, DynamicsAlgorithm, ClusteringAlgorithm
from music21 import scale
import numpy as np

from core import gen_random_seqs, load_seq_config

import os


def run():
    msa = load_seq_config(GLOBALS['TEST_VECTORS'] + '/two_gapped_seqs.json')
    print('MSA', msa)

    pitch_algo = PitchAlgorithm(PitchAlgorithm.WORD_DISTANCES, scale=scale.MajorScale().getPitches(), n_nucleotides=1)
    durations_algo = DurationsAlgorithm(DurationsAlgorithm.FREQUENCIES_DYNAMIC, window_size=10, window_duration=10, n_nucleotides=1)
    dynamics_algo = DynamicsAlgorithm(DynamicsAlgorithm.SHANNON_INDEX, window_size=10, gap_window_threshold=0.8, gap_column_threshold=0.9, criteria='local', levels=7)
    instruments_algo = ClusteringAlgorithm('kmeans')


    composer = Composer(instruments_algo, pitch_algo, durations_algo, dynamics_algo, input_type='array',
                        alignment=msa)

    sequence_names = np.array([str(i) for i in range(0, len(msa))])

    scores = composer.gen_numerical_vectors(k=2, piece_length=20)
    idx = 0
    for score in scores:
        fw = FileWriter(score, sequence_names)

        fname = 'demo_' + str(idx)

        print('fname', fname)
        fw.write(midi=fname, audio=fname, display=False)
        idx += 1

if __name__ == '__main__':
    run()