from ensemble import PitchAlgorithm, DurationsAlgorithm, DynamicsAlgorithm, ClusteringAlgorithm, Composer, FileWriter

from algorithms import *
import numpy as np

print('### Tests ###\n\n')
print('Test 1\n') 

msa = np.zeros((2,200), dtype="S1")

msa[0] = np.zeros(200, dtype="S1")

msa[0][0:50] = 'a'
msa[0][50:75] = 'c'
msa[0][75:150] = 'a'
msa[0][150:] = 'g'

print(msa[0])
msa[1] = np.array(msa[0], copy=True)
msa[1][75:100] = 'g'
msa[1][150:175] = 'a'

print(msa[1])
""""
aln = np.empty((len(msa), len(list(msa[0])[0])), dtype="S1")
for i in range(0, len(msa)):
    seq = list(msa[i])[0]
    for j in range(0, len(seq)):
        aln[i][j] = seq[j]

msa = aln
"""""

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