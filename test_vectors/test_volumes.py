import os

from algorithms import *
from config import GLOBALS
from ensemble import Composer, FileWriter

from core import load_seq_config

msa = load_seq_config(GLOBALS['TEST_VECTORS'] + '/two_gapped_seqs.json')
print msa
window_size = 20
window_duration = 20
n_nucleotides = 1

subdir = 'test_volume_demo'
os.environ['DIR'] = subdir

pitch_algo = PitchAlgorithm(PitchAlgorithm.WORD_DISTANCES, scale=scale.MajorScale().getPitches(), n_nucleotides=1)
durations_algo = DurationsAlgorithm(DurationsAlgorithm.FREQUENCIES_DYNAMIC, window_size=window_size, window_duration=window_duration, n_nucleotides=1)
dynamics_algo = DynamicsAlgorithm(DynamicsAlgorithm.SHANNON_INDEX, window_size=10, gap_window_threshold=0.3, gap_column_threshold=0.3, criteria='local', levels=7)
instruments_algo = ClusteringAlgorithm('kmeans')

composer = Composer(instruments_algo, pitch_algo, durations_algo, dynamics_algo, input_type='array',
                    alignment=msa)

sequence_names = np.array([str(i) for i in range(0, len(msa))])

i = 0

scores = composer.gen_numerical_vectors(k=2, piece_length=window_size*2)

print 'Parameters ' + str(window_size) + ' ' + str(window_duration) + ' ' + str(n_nucleotides)
print 'Number of scores', len(scores)
for score in scores:
    fw = FileWriter(score, sequence_names)
    fname = 'demo_' + str(i)

    print('fname', fname)
    fw.write(midi=fname, audio=fname, display=False, subdir=subdir)
    i += 1
    print('After writing', fname)
