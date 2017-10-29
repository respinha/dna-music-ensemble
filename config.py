from music21.instrument import StringInstrument, WoodwindInstrument, PitchedPercussion, BrassInstrument
import os

CURR_DIR = os.path.dirname(os.path.abspath(__file__))
SEQ_DIR = CURR_DIR + "/source_sequences"
OUTPUT_FILES = CURR_DIR + '/output_files'
MIN_TEMPO = 0.0625  # 64th

aln = lambda n: SEQ_DIR + '/mafft_' + str(n) + '.fasta'

GLOBALS = {'MEME_URL' : 'http://meme-suite.org/opal2/services/MEME_4.11.2',
           'SUPPORTED ALGORITHMS' : ['clustal', 'mafft', 'muscle'],
           'MAPPINGS' : {0: 'NO_SYNC', 1: 'SYNC_window_duration',2:'SYNC_BLOCL_DURATION_DISCRETE',3:'DURATION_WITH_DISTANCES'},
           'ALPHABET' : ['a', 'c', 'g', 't', '-'],
           'MUSE_SCORE' :   '/usr/bin/mscore',
           'FAMILIES' : {0: StringInstrument, 1: WoodwindInstrument,
                                  2: BrassInstrument, 3: PitchedPercussion},
           'SCORES' :   OUTPUT_FILES + '/scores',
           'MIDI' :   OUTPUT_FILES + '/midi',
           'AUDIO'  :   OUTPUT_FILES + '/audio',
           'HIST_DURATIONS' :   OUTPUT_FILES + '/stats/durations',
            'HIST_NOTES' :   OUTPUT_FILES + '/stats/notes',
           'ALIGNMENT_PARAMS' : ['fasta_file', 'seq_vector', 'n_seq', 'algorithm'],
           'TEST_VECTORS': CURR_DIR + '/test_vectors',
           'REGIONS_DIR': OUTPUT_FILES + '/regions'
           }
