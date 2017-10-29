import numpy as np
import matplotlib.pyplot as plt

from collections import defaultdict

from Bio import AlignIO
from Bio import SeqIO, Phylo

from Bio.Align import MultipleSeqAlignment

from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor

from music21 import instrument, note, stream, scale
from music21.instrument import StringInstrument, WoodwindInstrument, BrassInstrument, Percussion, PitchedPercussion
from music21 import midi

import itertools

from algorithms import *

import sys
import os
import shutil

from config import GLOBALS, OUTPUT_FILES, SEQ_DIR


class Composer(object):

    alignment = None

    def __init__(self, cluster_algorithm, pitch_algorithm,
                 duration_algorithm, dynamics_algorithm,
                 input_type='alignment', **kwargs):

        assert isinstance(cluster_algorithm, ClusteringAlgorithm) \
               and isinstance(pitch_algorithm, PitchAlgorithm) \
               and isinstance(duration_algorithm, DurationsAlgorithm)\
               and isinstance(dynamics_algorithm, DynamicsAlgorithm)\
               and input_type

        if input_type in ['alignment', 'array']:

            assert 'alignment' in kwargs.keys()
            self.alignment = kwargs['alignment']

        elif input_type == 'sequences':

            assert 'fasta' in kwargs.keys()

            seq_file = kwargs['fasta_file']

            if 'algorithm' not in kwargs.keys():
                print('Missing alignment algorithm')

            aln_algorithm = kwargs['algorithm']

            seq_vector, n_seq = None, None

            if 'seq_vector' not in kwargs.keys():

                if 'n_seq' not in kwargs.keys():
                    print('Missing explicit vector of sequences or maximum number of sequences')
                    sys.exit(1)

                n_seq = kwargs['n_seq']

            else:
                seq_vector = kwargs['seq_vector']

            from core import gen_alignment
            self.alignment = gen_alignment(seq_file, seq_vector=seq_vector, n_sequences=n_seq, algorithm=aln_algorithm)

        self.clustering_algorithm = cluster_algorithm
        self.durations_algorithm = duration_algorithm
        self.pitch_algorithm = pitch_algorithm
        self.dynamics_algorithm = dynamics_algorithm

    def assign_instruments(self):

        # assert isinstance(self.alignment, MultipleSeqAlignment)
        return [instrument.Timpani(), instrument.Glockenspiel(), instrument.Vibraphone(), instrument.Marimba()]

    def gen_numerical_vectors(self, k=2, piece_length=5):

        msa = AlignIO.read(self.alignment, 'clustal') if not isinstance(self.alignment, np.ndarray) else self.alignment

        window_sizes = np.zeros(3)
        if 'windows_size' in self.dynamics_algorithm.keys():
            window_sizes[0] = self.clustering_algorithm['windows_size']

        if 'windows_size' in self.pitch_algorithm.keys():
            window_sizes[1] = self.pitch_algorithm['windows_size']

        if 'windows_size' in self.durations_algorithm.keys():
            window_sizes[2] = self.durations_algorithm['windows_size']

        iterator = iter(window_sizes)
        window_size = next(iterator)
        window_sizes_are_valid = np.all(window_size == rest for rest in iterator if rest != 0)

        assert window_sizes_are_valid, 'Window sizes cannot differ in algorithm mappings'

        from core import gen_song

        # TODO: insert clustering/instrument assigning algorithm
        instruments = self.assign_instruments()[0:len(msa)]

        songs = gen_song(self.pitch_algorithm, self.durations_algorithm, self.dynamics_algorithm, msa, instruments, k, piece_length=piece_length)

        assert isinstance(songs, list)
        return songs

        #  dynamics_vector = gen_dynamics_vector(msa, self.dynamics_algorithm)
        # add_dynamics_to_score(dynamics_vector['vol'], score)

    def play(self):
        pass

    def shuffle_instruments(self):
        pass


class FileWriter(object):

    def __init__(self, score, records):

        self.score = score
        self.records = records

        assert isinstance(score, stream.Score) and isinstance(records, np.ndarray)

        print(score.parts, records)
        assert len(score.parts) == len(records), str(len(score.parts)) + ' ' + str(len(records))

    def write(self, name='alignment', display=False, stats=None, **kwargs):

        if len(kwargs.keys()) == 0:
            print('No output type or path specified')
            sys.exit(1)

        if not os.path.isdir(name):
            #shutil.rmtree(name)
            os.mkdir(name)

        subdir = None
        if 'subdir' in kwargs.keys():
            subdir = kwargs['subdir']

        # parsing arguments
        for key, value in kwargs.items():

            print(type(value))

            if key == 'midi':

                print('Writing midi...')

                if not value.endswith('.mid'):
                    value += '.mid'

                if not os.path.isdir(GLOBALS['MIDI']):
                    os.mkdir(GLOBALS['MIDI'])

                output_midi = GLOBALS['MIDI'] + '/'
                if subdir:
                    #    shutil.rmtree(output_midi + subdir)
                    output_midi += subdir + '/'
                    if not os.path.isdir(output_midi):
                        os.mkdir(output_midi)

                output_midi += value

                f = midi.translate.streamToMidiFile(self.score)

                print('Output MIDI', output_midi)

                f.open(output_midi, attrib='wb')
                f.write()
                f.close()

                import subprocess

                # currently audio can only be generated from MIDI input
                if 'audio' in kwargs.keys():
                    audio_name = output_midi.split('.mid')[0] + '.ogg'
                    print('AUDIO NAME', audio_name)
                    subprocess.call(["timidity", output_midi, "-Ow", "-o", audio_name])

            elif key == 'score':

                print('Writing score...')

                # lily already inserts extension in filename
                value = value[:-4] if value.endswith('.pdf') else value

                path = GLOBALS['SCORES'] + '/'
                if subdir:
                    #    shutil.rmtree(path + subdir)
                    path += subdir + '/'
                    if os.path.isdir(path):
                        os.mkdir(path)

                path += value
                self.score.write('lily.pdf', fp=path)

                # deleting tmp file created by lily without .pdf extension
                to_unlink = GLOBALS['SCORES'] + '/' + value
                if os.path.isfile(to_unlink):
                    os.unlink(to_unlink)

            if display:
                print('Displaying score...')
                self.score.show(app=GLOBALS['MUSE_SCORE'])

            if stats:

                assert isinstance(stats, str), 'Invalid path for stats file'

                print('Generating pitch and duration histograms...')

                hist_durations_dir = GLOBALS['HIST_DURATIONS']
                if not os.path.isdir(hist_durations_dir):
                    os.makedirs(hist_durations_dir)

                hist_notes_dir = GLOBALS['HIST_NOTES']
                if not os.path.isdir(hist_notes_dir):
                    os.makedirs(hist_notes_dir)

                i = 0
                for part in self.score.parts:

                    durations_notes = np.array([(float(x.duration.quarterLength), x.name) for x in part
                                                if isinstance(x, note.Note)],
                                               dtype=[('durations', np.float,), ('notes', 'S5')])

                    durations = durations_notes['durations']  # first column

                    print('Durations', durations)
                    print('Durations set', np.unique(durations))
                    notes = durations_notes['notes']

                    plt.hist(durations, color='b')
                    plt.savefig(hist_durations_dir + '/' + stats + '_' + self.records[i] + '.png')
                    plt.close()

                    from collections import Counter

                    counts = Counter(notes)

                    note_names = counts.keys()
                    note_values = counts.values()

                    width = 1.0
                    idx = np.arange(len(note_names))
                    plt.bar(idx, note_values, width)
                    plt.xticks(idx + width * 0.5, note_names)

                    plt.savefig(hist_notes_dir + '/' + stats + '_' + self.records[i] + '.png')
                    plt.close()

                    i += 1

# TODO:
# check if test file with results already exist  2 mins
# if not:
#   test BioPython
#   test music21
#   check if lily is in music21 environment
#   check if muse score is in music21 environment
#   check if pyqt is installed
#   check if scipy is installed
# else:
#   check file content for success or failure
def config_environment():
    pass

def run():

    #### config_environment() ####

    # test vectors
    # hand-made file
    # piece-length = 15
    # window = 5

    msa = np.array(
        [
            ['AAAAAAAAAAAAAAAAAAAACCCCCCCCCCCCCCCACGCTGTCAA'],
            ['AAAAAAAAAAAAAAAAAAAACCCCCCCCCCCCCCCAGCATCTCAG'],
            ['AAACGAAACGAAACGAAACGCCCCCCCCCCCCCCCCAGCAGTCAA'],
            ['AAATCAAATCAAATCAAATCCCCCCCCCCCCCCCCTTAGCGCTTA']
        ]
    )

    aln = np.empty((len(msa), len(list(msa[0])[0])), dtype="S1")

    for i in range(0, len(msa)):

        seq = list(msa[i])[0]
        for j in range(0, len(seq)):
            aln[i][j] = seq[j]

    msa = aln

    pitch_algo = PitchAlgorithm(PitchAlgorithm.WORD_DISTANCES, scale=scale.MajorScale().getPitches(), n_nucleotides=1)
    durations_algo = DurationsAlgorithm(DurationsAlgorithm.FREQUENCIES_DYNAMIC, window_size=2, window_duration=10, n_nucleotides=1)
    dynamics_algo = DynamicsAlgorithm(DynamicsAlgorithm.SHANNON_INDEX, window_size=2, gap_window_threshold=0.5, gap_column_threshold=0.7, criteria='local', levels=7)
    instruments_algo = ClusteringAlgorithm('kmeans')

    composer = Composer(instruments_algo, pitch_algo, durations_algo, dynamics_algo, input_type='array',
                        alignment=msa)

    sequence_names = np.array([str(i) for i in range(0, len(msa))])

    i = 0

    scores = composer.gen_numerical_vectors(k=2)
    for score in scores:
        fw = FileWriter(score, sequence_names)
        fname = 'hardcoded_test_' + str(i)

        print('fname', fname)
        fw.write(midi=fname, audio=fname, display=False)
        i += 1

    """aln_file = SEQ_DIR + '/clustal3.aln'

    from music21 import scale

    pitch_algo = PitchAlgorithm(PitchAlgorithm.WORD_DISTANCES, scale=scale.MajorScale().getPitches(), n_nucleotides=2)
    durations_algo = DurationsAlgorithm(DurationsAlgorithm.FREQUENCIES_DYNAMIC, window_size=1000, window_duration=500, n_nucleotides=2)
    dynamics_algo = DynamicsAlgorithm(DynamicsAlgorithm.SHANNON_INDEX, window_size=2500, gap_window_threshold=0.5, gap_column_threshold=0.7, criteria='local', levels=7)

    instruments_algo = ClusteringAlgorithm('kmeans')

    composer = Composer(instruments_algo, pitch_algo, durations_algo, dynamics_algo, input_type='alignment', alignment=aln_file)

    msa = AlignIO.read(aln_file, 'clustal')
    sequences = np.chararray(len(msa), itemsize=25)

    i = 0
    for record in msa:
        sequences[i] = record.description
        i += 1

    i = 0

    for k in range(2, 7):

        scores = composer.gen_numerical_vectors(k=k)

        for score in scores:

            fw = FileWriter(score, sequences)
            fname = 'test_fast_' + str(i) + "_" + str(k)

            print 'fname', fname
            fw.write(midi=fname, audio=fname, display=False)
            i += 1"""

    print('Done')

if __name__ == "__main__":
    
    import multiprocessing
    import time

    p = multiprocessing.Process(target=run, name="Run", args=())
    p.start()

    print('Started!')
    time.sleep(20 * 60)

    if p.is_alive():
        print('Early kill')
        p.terminate()
        p.join()