import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from Bio import Cluster

from collections import defaultdict

from Bio.SeqFeature import  SeqFeature
from Bio import AlignIO
from Bio import SeqIO, Phylo

from Bio.Align import MultipleSeqAlignment
from Bio.Phylo.PhyloXML import Phylogeny
from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor
from Bio.Phylo import Consensus
from Bio.Phylo.Applications import PhymlCommandline
from __builtin__ import object

from music21 import instrument, note, harmony, interval, key, duration
from music21 import stream, defaults

import itertools

from music21.instrument import StringInstrument, WoodwindInstrument, BrassInstrument, Percussion
from xcb.shape import shapeExtension

from api import GLOBALS

import sys


# abstract algorithms class
# cannot be instantiated by itself
class Algorithm(object):

    def __init__(self, algorithm):
        self.name = algorithm

    def __repr__(self):
        return self.name


# pitch algorithms
# supported types: [ dummy, discrete_sync, patternized ]
class DurationsAlgorithm(Algorithm):

    assures_sync = False

    def __init__(self, name='frequencies', sync=True):
        super.__init__(self)

        assert name is not None
        assert sync is not None

        self.name = name
        self.assures_sync = sync


# pitch algorithms
# supported types: [ word_distances, static_assign ]
class PitchAlgorithm(Algorithm):

    def __init__(self, name='word_distances', pitch_scale=None):
        super.__init__(self)

        assert name is not None
        self.name = name


# dynamics algorithms
# supported types: [ 'gaps' ]
class DynamicsAlgorithm(Algorithm):

    def __init__(self, name='gaps'):

        super.__init__(self)

        assert name is not None
        self.name = name


# species clustering
class ClusteringAlgorithm(Algorithm):

    def __init__(self, algorithm):
        super(ClusteringAlgorithm, self).__init__(algorithm)

algorithms = {'instruments':['kmeans','hierarchical'],'pitch':['distances'], 
                   'dynamics':['gaps'],'duration':['frequencies']}


class Ensemble(object):

    def __init__(self, mapper, alignment_file=True, **kwargs):

        assert isinstance(mapper, Mapper)
        assert mapper.ready_for_implementation()

        if alignment_file:

            assert 'alignment' in kwargs.keys()
            self.alignment = kwargs['alignment']
        else:

            if 'fasta_file' not in kwargs.keys():
                print('Missing input file for alignment')
                sys.exit(0)

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

            from api import gen_alignment
            self.alignment = gen_alignment(seq_file, seq_vector=seq_vector, n_sequences=n_seq, algorithm=aln_algorithm)

        self.mapper = mapper

        """for x, y in mapper.items():

            assert y in algorithms[x]

            if x == 'instruments':

                self.instrument_vector = ClusteringAlgorithm(y)

            elif x == 'pitch':

                self.pitch_vector = PitchAlgorithm(y)

            elif x == 'duration':

                self.durations_vector = DurationsAlgorithm(y)

            elif x == 'dynamics':

                self.dynamics_vector = DynamicsAlgorithm(y)"""

    def assign_instruments(self, msa, nclusters=0):

        assert isinstance(msa, MultipleSeqAlignment)

        algorithm = self.mapper['instruments']

        from api import get_clusters_from_alignment
        return get_clusters_from_alignment(msa, algorithm=algorithm, nclusters=nclusters)

    def gen_numerical_vectors(self, **kwargs):

        msa = AlignIO.read(self.alignment)

        window = 0
        if 'window' in kwargs.keys():

            window = kwargs['window']
            assert isinstance(window, int) or isinstance(window, float)

        window_duration = 0
        if 'window_duration' in kwargs.keys():

            window_duration = kwargs['window_duration']
            assert isinstance(window_duration, int) or isinstance(window_duration, float)

        from api import numeric_vectors, dynamics_algorithm

        instruments = None
        if 'nclusters' in kwargs.keys():

            nclusters = kwargs['nclusters']
            self.assign_instruments(msa, nclusters=nclusters)
        else:
            self.assign_instruments(msa)

        if 'dynamics_window' not in kwargs.keys():
            dynamics_window = None
        else:
            dynamics_window = kwargs['dynamics_window']

        return instruments, dynamics_algorithm(msa, algorithm='entropy', window=dynamics_window),\
               numeric_vectors(msa, block=window,
                               block_duration=window_duration,
                               duration_mapping=self.mapper['duration'],
                               )
    def play(self):
        pass

    def shuffle_instruments(self):
        pass

    def write_file(self, filename=None):
        fw = FileWriter()
        pass


class Mapper(dict):

    def __init__(self, **kwargs):

        super(Mapper, self).__init__(**kwargs)

        for x,y in kwargs.items():

            if x not in algorithms:
                print(x + ' is not a valid algorithm type')
            else:
                self[x] = y

    def __getattr__(self, attr):
        return self[attr]

    def __setattr__(self, key, value):

        assert isinstance(value, Algorithm)
        assert key in algorithms

        self[key] = value

    def ready_for_implementation(self):

        from collections import Counter
        return Counter(algorithms) == Counter(self.keys())


class FileWriter(object):

    def __init__(self, ensemble):

        assert ensemble is not None and isinstance(ensemble, Ensemble)
        self.ensemble = ensemble

def run():

    algorithms_ = {'instruments':'kmeans','pitch':'','duration':'','dynamics':''}

    mapper = Mapper(algorithms_)
    ensemble = Ensemble(mapper)

    fw = FileWriter(ensemble)

    # app.run()
    pass
""""
class Species(object):

    def __init__(self, sequence, mapping=default_mapping):

        if mapping == default_mapping:
            self._mapping = MusicDNAMapping()
            self.track = self._mapping.process_sequence(sequence)

        else:
            raise MappingException('No valid DNA-music mapping was introduced')


class Orchestra(object):

    instrument_families = {0: StringInstrument, 1: WoodwindInstrument,
                           2: BrassInstrument, 3: Percussion}

    instruments = []

    def __init__(self, alignment, mapping=default_mapping):

        assert alignment is not None and isinstance(alignment, file) and os.path.isfile(alignment)
        assert mapping is not None and mapping in GLOBALS['VALID_MAPPING']

        self._alignment = alignment

        clusters = cluster_alignment(alignment)

        for c in clusters:
            self.gen_instrument(c)

        alignment = alignment.split('.')[0] + '.phy'
        alignment = AlignIO.read(open(alignment), 'rU')

        summary = AlignInfo.SummaryInfo(alignment)
        expect_freq = {'A': .25, 'G': .25, 'T': .25, 'C': .25}

        from Bio.Alphabet import IUPAC
        from Bio.SubsMat import FreqTable

        e_freq_table = FreqTable.FreqTable(expect_freq, FreqTable.FREQ, IUPAC.unambiguous_dna)

        # entropy estimation
        while i < len(alignment[0]):

            # perform musical mapping
            h = summary.information_content(i, i + 20, e_freq_table=e_freq_table, chars_to_ignore=['N', '-'])
            # module musical mapping with entropy value

    def gen_instrument(self, family_number):

        family = self.instrument_families[family_number]
        subclasses = family.__subclasses__()

        # totally random for now
        cl = os.urandom(len(subclasses)-1)
        cl = subclasses[cl]

        self.instruments.append(cl)


    def create_midi(self, mapping=default_mapping):

        for track in self.tracks:
            print 'track'

            self.composition.set_author("respinha", "ruiespinharibeiro@gmail.com")
            self.composition.set_title("Really Dummy Composition", "Virtual Mitochondria Orchestra")
            self.composition.add_track(track)

            print self.composition.selected_tracks

        #midi_file_out.write_Composition('orchestra.mid', composition=self.composition)
        #midi_file_out.write_Composition('composition', LilyPond.from_Composition(self.composition))


class MusicDNAMapping(object):

    def __init__(self, mapping=default_mapping):
        self._notes = dict()

        self.mapping = mapping
        if mapping == default_mapping:

            self._notes['C'] = note.Note("C4")
            self._notes['A'] = note.Note("C4")
            self._notes['G'] = note.Note("C4")
            self._notes['T'] = note.Note("C", 5)
            self._notes['-'] = note.Rest(type='whole')

            #self._bar = Bar()    # default tempo, 4/4
            self._bar = []

    def process_sequence(self, sequence):
        
        length = len(sequence)
        
        #bars = np.array(self._bar for _ in range(length))
        bars = [self._bar for _ in range(0, length)]

        for i in range(0, len(sequence)):

            char = str(sequence[i])
            if char in self._notes.keys():
                bars[i%4].place_notes(self._notes[char], 4)

        #t = Track(Violin())
        #for b in bars: t.add_bar(b)

        #return t

        return None

    def frequency_vectors(self, sequence):

        frequency_vect = dict()

        step = 100
        for inc in range(0, len(sequence), step):

            subseq = sequence[inc-step:inc]
            #print subseq.seq

            series = pd.Series(subseq)
            counts = dict(series.value_counts(normalize=True))

            #print counts
            for key, value in counts.iteritems():
                if key not in frequency_vect:
                    frequency_vect[key] = [value]
                else:
                    frequency_vect[key].append(value)

                print frequency_vect

        for key, val in frequency_vect.iteritems():
            frequency_vect[key] = np.array(val)
        return frequency_vect

    def distance_vectors(self, sequence):

        vectors = dict() # keys are nucleotides; values are np arrays

        last_occurrence = dict()

        if self.mapping == 'Distances':

            length = len(sequence)
            for i in range(0, length):

                letter = sequence[i]

                if letter not in vectors.keys():

                    #vectors[letter] = np.zeros(shape=(length))
                    vectors[letter] = []
                else:
                    diff = i - last_occurrence[letter]
                    vectors[letter].append(diff)

                last_occurrence[letter] = i

            for x, y in vectors.iteritems():
                diff = length - last_occurrence[x]
                vectors[x].append(diff)

            return vectors

        else:
            raise MappingException('DNA Mapping ' + self.mapping + ' does not allow numeric transformation')


class MappingException(Exception):

    def __init__(self, message):
        super(MappingException, self).__init__(message)

    def __repr__(self):
        return self.message


if __name__ == '__main__':

    print 'Clustering from 4 to 7 using default algorithm (Clustal)'

    for i in range(7,11):
        cluster_alignment(aln(i), depth=2)

    #n = int(argv[2])

    #size = int(argv[1])
    #n = int(argv[2])

    #alignment_path = 'source_sequences/clustal_' + str(size) + ".aln"

    #orchestra = Orchestra(alignment, tree)
    #orchestra.create_midi()"""""