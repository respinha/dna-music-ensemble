from ensemble import Composer, FileWriter
from algorithms import PitchAlgorithm, DurationsAlgorithm, DynamicsAlgorithm, ClusteringAlgorithm
from music21 import scale
import numpy as np

print('### Tests ###\n\n')
print('Test 1\n')


def run():
    msa = np.zeros((2, 200), dtype="S1")

    #msa[0] = np.zeros(200, dtype="S1")
    msa[0] = np.random.choice(['a','g','c','t','-'], 200)
    print(msa[0])

    msa[1] = np.random.choice(['a','g','c','t','-'], 200)
    print(msa[1])

    error_log = open('test_vectors/error_log_test2.txt', 'wr')

    for window_size in range(5, 30, 5):
        for window_duration in range(5, 60, 5):
            for n_nucleotides in range(1, 2):
                try:
                    subdir = 'demo2_' + str(window_size) + '_' + str(window_duration) + '_' + str(n_nucleotides)

                    pitch_algo = PitchAlgorithm(PitchAlgorithm.WORD_DISTANCES, scale=scale.MinorScale().getPitches(),
                                                n_nucleotides=1)
                    durations_algo = DurationsAlgorithm(DurationsAlgorithm.FREQUENCIES_DYNAMIC, window_size=window_size,
                                                        window_duration=window_duration, n_nucleotides=1)
                    dynamics_algo = DynamicsAlgorithm(DynamicsAlgorithm.SHANNON_INDEX, window_size=window_size,
                                                      gap_window_threshold=0.5, gap_column_threshold=0.7,
                                                      criteria='local', levels=7)
                    instruments_algo = ClusteringAlgorithm('kmeans')

                    composer = Composer(instruments_algo, pitch_algo, durations_algo, dynamics_algo, input_type='array',
                                        alignment=msa)

                    sequence_names = np.array([str(i) for i in range(0, len(msa))])

                    i = 0

                    scores = composer.gen_numerical_vectors(k=2, piece_length=window_size * 2)

                    print 'Parameters ' + str(window_size) + ' ' + str(window_duration) + ' ' + str(n_nucleotides)
                    print 'Number of scores', len(scores)
                    for score in scores:
                        fw = FileWriter(score, sequence_names)
                        fname = 'demo_' + str(i)

                        print('fname', fname)
                        fw.write(midi=fname, audio=fname, display=False, subdir=subdir)
                        i += 1
                        print('After writing', fname)

                except Exception as e:
                    print 'Error ocurred! Please check log!'
                    error_log.write('Error occurred with windows_size = ' + str(window_size) + \
                                    ' & window_duration = ' + str(window_duration) + ' & n_nucleotides = ' + str(
                        n_nucleotides) + e.message + '\n')


if __name__ == '__main__':
    run()