from __future__ import division

import subprocess

import time
from Bio import AlignIO
from Bio import Phylo
from Bio import SeqIO
from Bio.Align import AlignInfo, MultipleSeqAlignment

from scipy.stats import itemfreq

import numpy as np

import datasketch as dk
import pandas as pd

import os
import sys
import random

from matplotlib import pyplot as plt

from music21 import note, scale, instrument, stream, duration, tempo
from music21.instrument import StringInstrument, WoodwindInstrument, BrassInstrument, PitchedPercussion, Instrument
from music21 import Music21Object

from algorithms import *
from config import GLOBALS, MIN_TEMPO
from core.bio import msa_to_phylip


def gen_dynamics_vector(msa, dynamics_algorithm):
    # criteria: local, avg ou median entropy
    assert isinstance(dynamics_algorithm, DynamicsAlgorithm)
    assert 'window_size' in dynamics_algorithm.keys(), 'Empty window for dynamics algorithm'
    assert dynamics_algorithm['algorithm'] == DynamicsAlgorithm.SHANNON_INDEX

    window = dynamics_algorithm['window_size']

    if 'gap_column_threshold' not in dynamics_algorithm.keys():
        gap_column_threshold = 0.7
    else:
        gap_column_threshold = dynamics_algorithm['gap_column_threshold']

    if 'gap_window_threshold' not in dynamics_algorithm.keys():
        gap_window_threshold = 0.7
    else:
        gap_window_threshold = dynamics_algorithm['gap_window_threshold']

    if 'criteria' not in dynamics_algorithm.keys():
        criteria = 'local'
    else:
        criteria = dynamics_algorithm['criteria']

    if 'levels' not in dynamics_algorithm.keys():
        levels = 5
    else:
        levels = dynamics_algorithm['levels']

    aln_len = len(msa[0])

    from math import ceil

    n_windows = np.int(ceil(float(aln_len) / window))

    gaps_below_thresh = np.zeros(n_windows, dtype=np.bool)

    window_idx = 0

    # first iterating through MSA to identify
    # which windows have a percentage of gaps
    # below the established threshold
    for i in range(0, aln_len, window):

        boundary = window \
            if i + window <= aln_len \
            else aln_len - i

        local_is_below_threshold = np.zeros((boundary,), dtype=np.bool)

        for j in range(i, i + boundary):

            column = np.array([c for c in msa[:, j]])

            n_gaps = np.count_nonzero(column == '-')

            if float(n_gaps) / len(column) < gap_column_threshold:
                local_is_below_threshold[j - i] = True

        n_ungapped_regions = len(np.where(local_is_below_threshold)[0])
        print 'N ungapped regions', n_ungapped_regions, gap_window_threshold * boundary
        if n_ungapped_regions >= gap_window_threshold * boundary:
            gaps_below_thresh[window_idx] = True
        window_idx += 1

    #print('Gaps below thresh', gaps_below_thresh)

    from scipy.stats import entropy

    dynamics_vector = np.zeros((n_windows,),
                               dtype=[('entropy', np.float), ('vol', np.float)])
    entropies_idx = 0

    for i in range(0, aln_len, window):
        boundary = window \
            if i + window <= aln_len \
            else aln_len - i

        # if this window has a percentage of gaps
        # above the considered threshold
        if not gaps_below_thresh[entropies_idx]:
            dynamics_vector['entropy'][entropies_idx] = -1
            entropies_idx += 1
            continue

        local_entropy = np.zeros((boundary,))

        for j in range(i, i + boundary):

            column = msa[:, j]
            column_symbols = column[np.where(column != '-')]

            if len(column_symbols) < (1 - gap_column_threshold) * len(column):
                # apply thresholding
                #pass
                local_entropy[j - i] = 0
            else:
                # Shannon entropy of all column symbols
                counts = itemfreq(column_symbols)
                counts = np.array([float(count) for count in counts[:, 1]])

                counts /= np.sum(counts)

                local_entropy[j - i] = entropy(counts, base=2)

        if criteria == 'local':
            dynamics_vector['entropy'][entropies_idx] = np.sum(local_entropy)
        elif criteria == 'average':
            dynamics_vector['entropy'][entropies_idx] = np.average(local_entropy)
        elif criteria == 'median':
            dynamics_vector['entropy'][entropies_idx] = np.median(local_entropy)
        else:
            print 'Unsupported criteria ' + str(criteria) + ' for entropy aggregation'
            sys.exit(1)

        entropies_idx += 1

    max_vol = 0.95
    min_vol = 0.30

    entropies = dynamics_vector['entropy']

    split_info = np.array_split(np.sort(np.unique(entropies)), levels)  # splitting info into classes
    volumes = np.linspace(min_vol, max_vol, num=levels)  # vector with all possible volumes

    for i in range(0, int(n_windows)):

        if entropies[i] == -1:
            dynamics_vector['vol'][i] = -1

        for j in range(0, len(split_info)):
            if entropies[i] <= split_info[j][-1]:
                dynamics_vector['vol'][i] = volumes[j]
                break

    return dynamics_vector


def add_dynamics_to_score(dynamics_vector, score, window_size, instruments, max_rest_tempo=3):

    # consistency check
    assert isinstance(score, stream.Score) and isinstance(dynamics_vector, np.ndarray)
    assert len(instruments) == len(score.parts)
    assert all(isinstance(i, instrument.Instrument) for i in instruments)

    score_tempo = score.getElementsByClass(tempo.MetronomeMark)[0]

    # creating filtered score
    # starting with empty parts
    final_score = stream.Score()
    length = np.inf

    for p in range(0, len(score.parts)):

        part_len = len(score.parts[p]) - 1 # excluding instrument elem!!
        if part_len < length:
            length = part_len

        part = stream.Part()
        part.insert(0, score_tempo)
        part.insert(0, instruments[p])

        final_score.append(part)

    # iterating through a music score in chunks of 'window_size' length
    vol_idx = 0

    for i in range(0, length, window_size):
        boundary = window_size if i + window_size <= length else length - i

        if dynamics_vector[vol_idx] <= 0: # silence
            r = note.Rest()

            for part in final_score.parts:
                part.append(r)
                part[-1].seconds = max_rest_tempo

        else:
            # iterating over parts
            for j in range(0, len(final_score.parts)):

                # scores have same number of parts
                final_part = final_score.parts[j]
                part = score.parts[j]

                for k in range(i, i + boundary):

                    n = part[k]
                    if isinstance(n, note.GeneralNote):  # if Note or Rest

                        final_part.append(n)
                        if isinstance(n, note.Note):  # if Note

                            final_part[-1].volume = dynamics_vector[vol_idx]

        vol_idx += 1

    print 'SECS', final_score[0].seconds

    assert len(final_score.parts) == len(score.parts)
    return final_score


# returns a tuple containing:
#   - a word distance vector per word (A,C,G,T)
#   - a label array with the assigned duration of each nucleotide
def gen_pitch_duration_vectors(sequence, pitch_algorithm, durations_algorithm):
    assert isinstance(pitch_algorithm, PitchAlgorithm) and isinstance(durations_algorithm, DurationsAlgorithm)
    assert sequence is not None

    # step is the number of columns/characters that are mapped
    length = len(sequence)

    # auxiliary structures
    distance_vectors = dict()  # keys are nucleotides; values are np arrays
    last_occurrence = dict()  # aux dict used for word distances

    step = 1
    window = 1500
    window_duration = 8

    for key in durations_algorithm.keys():
        if key == 'n_nucleotides':
            step = durations_algorithm['n_nucleotides']
            assert isinstance(step, int) and step > 0

        elif key == 'window_size':
            window = durations_algorithm['window_size']
            assert isinstance(window, int) and window > 0

        elif key == 'window_duration':
            window_duration = durations_algorithm['window_duration']
            assert (isinstance(window_duration, float) or isinstance(window_duration, int)) and window > 0

    if 'n_nucleotides' in pitch_algorithm.keys():
        assert step == pitch_algorithm['n_nucleotides']

    d_algorithm = durations_algorithm['algorithm']
    p_algorithm = pitch_algorithm['algorithm']

    assert 'scale' in pitch_algorithm.keys()

    if p_algorithm == PitchAlgorithm.WORD_FREQ:
        scale = p_algorithm['scale']
        pitch_freq_vectors = dict()

        assert len(scale) > 0

    # splitting by window value
    split_seq_len = int(length / window) if length % window == 0 else int(length / window) + 1
    split_sequence = np.zeros((split_seq_len, window), dtype="S1")

    for i in range(0, split_seq_len):
        subseq = sequence[i * window: i * window + window]
        split_sequence[i][: len(subseq)] = subseq
        split_sequence[len(subseq):] = ''

    sequence = split_sequence

    n_nucleotide_split_len = int(window / step) if window % step == 0 else int(window / step) + 1
    split_sequence = np.zeros((split_seq_len, n_nucleotide_split_len), dtype="S" + str(step))

    nucleotides_idx = 0
    for i in range(0, len(sequence)):

        subseq = sequence[i]
        split_subseq = split_sequence[i]

        for j in range(0, len(subseq), step):

            for k in range(j, j + step):
                split_subseq[nucleotides_idx] += subseq[k]

            nucleotides_idx += 1

        nucleotides_idx = 0

    # sequence shape - [ [window] [window] ....]
    # window shape - ['n-nucleotide' 'n-nucleotide' ....]
    sequence = split_sequence

    # print sequence
    assert len(sequence.shape) == 2

    offset = 0

    # 1d ndarray containing durations in quarter lengths
    durations = np.empty((sequence.shape[0] * sequence.shape[1],), dtype=np.float)

    for i in range(0, len(sequence)):

        # a window with n-nucleotide elements
        subset = sequence[i]
        # grouping frequencies of n-nucleotides  within a block
        counts = dict(itemfreq(subset))

        # This algorithm assigns durations dynamically to nucleotides
        # based on their relative frequency
        if d_algorithm == durations_algorithm.FREQUENCIES_DYNAMIC or p_algorithm == PitchAlgorithm.WORD_FREQ:

            counts_sum = 0

            # for j in range(i, bounverifiyingdary):
            for j in range(0, len(subset)):
                counts_sum += int(counts[subset[j]])

        total_time = 0.0

        for j in range(0, len(subset)):

            if subset[j] == '': continue

            local_count = int(counts[subset[j]])
            freq = float(local_count) / len(subset)

            # pitch algorithm
            if p_algorithm == PitchAlgorithm.WORD_DISTANCES or d_algorithm == DurationsAlgorithm.WORD_DISTANCES:

                # TODO: optimization
                #   fazer uma primeira passagem pela sequencia inteira
                #   para saber o numero exato de nucleotidos e reservar np.arrays??
                letter = subset[j]

                if letter not in distance_vectors.keys():
                    distance_vectors[letter] = []

                else:

                    diff = j + offset - last_occurrence[letter]
                    distance_vectors[letter].append(diff)

                last_occurrence[letter] = j + offset

            if p_algorithm == PitchAlgorithm.WORD_FREQ:

                # e.g. major scale - C-D-E-F-G-A-B-C
                # obtain scale from p_algorithm
                # assign fraction to each note of scale
                # assign fraction of a note to a frequency

                len_scale = len(scale)
                frac = 1.0 / len(scale)

                for k in range(0, len_scale):
                    val = frac * k
                    if freq > val:
                        continue
                    else:

                        if letter not in pitch_freq_vectors.keys():
                            pitch_freq_vectors[letter] = []

                        pitch_freq_vectors[letter].append(val)

            # durations algorithm
            # frequency-biased algorithms
            if d_algorithm == DurationsAlgorithm.FREQUENCIES_DYNAMIC:

                assert counts_sum is not None and counts_sum > 0, 'Inconsistent values for counts of characters in region'

                local_duration = float(window_duration) * float(local_count) / float(counts_sum)
                total_time += local_duration

                # assert local_duration > MIN_TEMPO, \
                #    'Higher tempo required for each subsequence; too short duration was calculated: ' + str(local_duration
                if local_duration < MIN_TEMPO:
                    local_duration = MIN_TEMPO

            # duration biased with discrete attributions
            elif d_algorithm == DurationsAlgorithm.FREQUENCIES_DISCRETE:

                duration_labels = np.array(['32nd', '16th', 'eighth', 'quarter', 'half', 'whole'])
                frequencies = np.linspace(0.0, 0.5, num=len(duration_labels) - 1)

                # dictionary mapping a discrete value for word frequencies to each duration label
                duration_labels = {frequencies[k]: duration_labels[k] for k in range(0, len(duration_labels) - 1)}
                duration_labels['whole'] = 1.0

                local_duration = '32nd'
                keys = duration_labels.keys()

                # finding correspondent label
                for k in keys:
                    if k > freq:
                        local_duration = duration_labels[k]
                        break

                # obtaining duration in quarter lengths
                local_duration = duration.Duration(local_duration).quarterLength
                total_time += local_duration

            elif d_algorithm == DurationsAlgorithm.WORD_DISTANCES:

                duration_labels = np.array(['64th', '32nd', '16th', 'eighth', 'quarter', 'half', 'whole'])

                if len(distance_vectors[letter]) > 0:
                    duration_label = duration_labels[distance_vectors[letter][-1] % len(duration_labels)]
                else:
                    duration_label = 'eighth'

                local_duration = duration.Duration(duration_label).quarterLength
                total_time += local_duration

            else:
                print 'Invalid mapping introduced: ', d_algorithm['algorithm']
                raise NotImplementedError

            # durations[j] = local_duration
            durations[j + offset] = local_duration

        if d_algorithm == DurationsAlgorithm.FREQUENCIES_DISCRETE or d_algorithm == DurationsAlgorithm.WORD_DISTANCES:

            if round(float(total_time), 5) != round(float(window_duration), 5):

                ratio = float(window_duration) / float(total_time)
                total_time = 0.0
                # for j in range(i, boundary):
                for j in range(0, len(subset)):
                    durations[j + offset] *= ratio

                    assert durations[j + offset] >= MIN_TEMPO, \
                        'Higher tempo required for each subsequence; too short duration was calculated: ' + str(
                            durations[j + offset])

                    total_time += durations[j + offset]

        offset += len(subset)

    if p_algorithm == PitchAlgorithm.WORD_DISTANCES:

        for x in distance_vectors.keys():  # iter(distance_vectors.items()):

            diff = offset + 1 - last_occurrence[x]
            distance_vectors[x].append(diff)

        d_vectors_len = 0

        for x in distance_vectors.keys():
            distance_vectors[x] = np.array(distance_vectors[x])
            d_vectors_len += len(distance_vectors[x])

            # assert d_vectors_len * step == length, "Lengths don't match: sequence length = " + str(length) + "; d_vectors length: " + str(d_vectors_len)
    else:

        return pitch_freq_vectors, durations

    # (distances, frequencies per N words)
    print 'PITCH DURATION VECTORS', distance_vectors, durations
    return distance_vectors, durations


def gen_stream(score, sequence, pitch_algorithm, durations_algorithm, assigned_instrument):
    assert isinstance(pitch_algorithm, PitchAlgorithm) and isinstance(durations_algorithm, DurationsAlgorithm)

    if 'window_size' in durations_algorithm.keys():
        assert len(sequence) > durations_algorithm['window_size'], \
            'Invalid piece and window size ' + str(len(sequence)) + ' ' + str(durations_algorithm['window_size'])

    dv, durations = gen_pitch_duration_vectors(sequence, pitch_algorithm, durations_algorithm)

    for x in dv.keys():
        dv[x] = iter(dv[x])

    durations = iter(durations)

    # TODO: integrar esta parte no algoritmo anterior para poupar iteracoes
    # scale_len = len(scale.MajorScale().getPitches())

    assert 'scale' in pitch_algorithm.keys()
    s = pitch_algorithm['scale']
    scale_len = len(s)

    assert isinstance(s, list) and len(s) > 1

    part = stream.Part()

    # assert isinstance(score_tempo, tempo.MetronomeMark) and score_tempo == score.getElementsByClass(tempo.MetronomeMark)[0]
    print 'Assigned instrument', assigned_instrument
    part.insert(0, assigned_instrument)

    print 'Assigning notes and durations from numeric vectors...'

    # for l in sequence:
    step = pitch_algorithm['n_nucleotides']
    window = durations_algorithm['window_size']

    # TODO: make method
    # splitting by window value
    length = len(sequence)
    split_seq_len = int(length / window) if length % window == 0 else int(length / window) + 1
    split_sequence = np.zeros((split_seq_len, window), dtype="S1")

    for i in range(0, split_seq_len):
        subseq = sequence[i * window: i * window + window]
        split_sequence[i][: len(subseq)] = subseq
        split_sequence[len(subseq):] = ''

    sequence = split_sequence

    # splitting in n-nucleotides
    n_nucleotide_split_len = window / step if window % step == 0 else window / step + 1
    # print 'split seq len', split_seq_len, 'n_nucleotide split len', n_nucleotide_split_len
    split_sequence = np.zeros((split_seq_len, np.int(n_nucleotide_split_len)), dtype="S" + str(step))

    nucleotides_idx = 0
    for i in range(0, len(sequence)):

        subseq = sequence[i]
        split_subseq = split_sequence[i]

        for j in range(0, len(subseq), step):

            for k in range(j, j + step):
                split_subseq[nucleotides_idx] += subseq[k]

            nucleotides_idx += 1

        nucleotides_idx = 0

    sequence = split_sequence

    #print 'Shapes', sequence.shape[0], sequence.shape[1]
    symbols = ''
    for i in range(0, len(sequence)):

        subseq = sequence[i]

        for symbol in subseq:

            if symbol == '': continue

            pitch = dv[symbol].next()
            d = durations.next()

            symbols += symbol
            # if symbol == '-': print 'hey', set(symbol)
            if list(set(symbol)) != ['-']:

                if pitch_algorithm['algorithm'] == PitchAlgorithm.WORD_DISTANCES:

                    n = pitch % scale_len
                    # n = s.getPitches()[n]
                    n = s[n]
                    n = note.Note(n)

                else:

                    frac = 1.0 / len(scale)

                    for i in range(0, len(scale)):
                        val = i * frac

                        if pitch == val:
                            n = s[i]
                            n = note.Note(n)

                n.duration = duration.Duration(d)

            else:
                n = note.Rest()
                n.duration = duration.Duration(d)

                n.addLyric(symbol)

                assert isinstance(n, Music21Object)

            part.append(n)

    print('Inserting part on score', len(part))
    score.insert(0, part)
    print('Done')


def gen_song(pitch_algorithm, durations_algorithm, dynamics_algorithm, alignment, instruments, k_shingles,
             piece_length=5000):
    ####### ALIGNMENT HANDLING ##############
    assert (alignment is not None), 'No MSA provided'

    assert isinstance(alignment, MultipleSeqAlignment) or \
           (isinstance(alignment, str) and os.path.isfile(alignment)) or \
           (isinstance(alignment, np.ndarray) and len(alignment.shape) == 2)

    if isinstance(alignment, str):
        print 'Reading alignment...'
        aln_file = AlignIO.read(open(alignment, 'rU'), 'clustal')
        aln_file = msa_to_phylip(aln_file)

        print 'Opening phylip file...'
        alignment = AlignIO.read(open(aln_file.split('.')[0] + ".phy"), 'phylip-relaxed')

    assert isinstance(piece_length, int) or isinstance(piece_length, float)  # and piece_length > 60

    # piece_length for now is only referring to number of musical elements
    # n_pieces = len(alignment[0]) / (step * piece_length)
    scores = []
    if not isinstance(alignment, np.ndarray):
        alignment = np.array([[y for y in x] for x in alignment])

    # k = np.random.choice(np.arange(3, 7), 1)[0] # random number between 3 and 6; used for k-shingling
    print 'K =', k_shingles

    from core.music.similarity import SimHandler

    split_len = int(alignment.shape[1] / piece_length) \
        if alignment.shape[1] % piece_length == 0 \
        else int(alignment.shape[1] / piece_length) + 1

    split_alignment = np.array_split(alignment, split_len, axis=1)

    sim = SimHandler(split_alignment, k=k_shingles)
    clusters = sim.cluster_by_similarites()

    print('Clusters', clusters)
    tempos = np.arange(45, 160, (160 - 45) / len(clusters))
    tempos_vector = sim.assign_tempos_by_clusters(clusters, tempos)

    assert len(tempos_vector) == len(clusters) == len(split_alignment)

    piece_idx = 0

    for p in range(0, alignment.shape[1], piece_length):

        if alignment.shape[1] - p < piece_length:
            piece_length = alignment.shape[1] - p

        score = stream.Score()

        score_tempo = tempo.MetronomeMark('tempo', tempos_vector[piece_idx])

        # print 'SCORE TEMPO', score_tempo._number
        score.insert(0, score_tempo)

        print 'Generating pitches and durations...'

        subsequence = alignment[:, p: p + piece_length]

        regions_file_path = 'regions_' + str(piece_idx) + '.txt'

        if not os.path.isdir(GLOBALS['REGIONS_DIR']):
            os.mkdir(GLOBALS['REGIONS_DIR'])

        if not 'DIR' in os.environ.keys(): os.environ['DIR'] = 'default'
        if not os.path.isdir(GLOBALS['REGIONS_DIR'] + '/' + os.environ['DIR']):
            os.mkdir(GLOBALS['REGIONS_DIR'] + '/' + os.environ['DIR'])

        regions_file_path = GLOBALS['REGIONS_DIR'] + '/' + os.environ['DIR'] + '/' + regions_file_path
        regions_file = open(regions_file_path, 'wr')

        for i in range(0, alignment.shape[0]):
            regions_file.write(''.join(subsequence[i]) + '\n')
            gen_stream(score, subsequence[i], pitch_algorithm, durations_algorithm, instruments[i])

        print 'Checking if parts have the same total duration...'

        # aligning part durations (score or midi cannot be produced with unequal duration parts)
        for part in score.parts:

            # obtaining highest duration from all parts
            # and aligning with it
            diff = score.highestTime - part.highestTime

            if diff > 0:

                while round(diff, 5) > 0:

                    # minimum = duration.Duration('2048th')
                    n = note.Rest()

                    if diff >= float(0.5):
                        n.duration = duration.Duration(0.5)
                    else:
                        if diff >= MIN_TEMPO:
                            n.duration = duration.Duration(diff)
                        else:
                            n.duration = duration.Duration(MIN_TEMPO)

                    assert n.duration.quarterLength >= MIN_TEMPO

                    part.append(n)
                    diff = score.highestTime - part.highestTime

        dynamics_vector = gen_dynamics_vector(subsequence, dynamics_algorithm)

        volumes = dynamics_vector['vol']
        print 'VOLUMES', dynamics_vector
        window_size = dynamics_algorithm['window_size']

        score = add_dynamics_to_score(volumes, score, window_size, instruments)

        print 'Dynamics to score'
        """for part in new_score:
            elems = ''
            for y in range(2, len(part)):
                elems += str(part[y].volume) + ' '
            print elems
        sys.exit(1)"""""
        scores.append(score)

        regions_file.write('\n\nTempo: ' + str(tempos_vector[piece_idx]))
        regions_file.close()

        piece_idx += 1

    """print similarities_df

    classes = np.arange(40, 145, (145 - 40) / n_classes, dtype=np.uint8)
    print classes

    jaccard_indices = np.array_split(similarities_df['jaccard'], n_classes)
    class_size = len(jaccard_indices[0])

    j = 0
    for i in range(0, len(classes)):
        similarities_df['tempo'][j:j + class_size] = classes[i]
        print 'L', len(similarities_df['tempo'][j:j + class_size])
        j += class_size

    similarities_df = similarities_df.sort_values('idx')"""

    # parte estatistica e output de ficheiros para @FileWriter
    # retornar score, utilizar dynamics_algorithm, adicionar volumes a score e analisar score
    return scores

"""""""""
if __name__ == "__main__":
    from config import CURR_DIR
    msa = AlignIO.read(CURR_DIR + '/output.fasta', 'clustal')
    msa = np.array([[y for y in x] for x in msa])

    window_size = 1500
    piece_length = 5000

    dynamics_algorithm = DynamicsAlgorithm(DynamicsAlgorithm.SHANNON_INDEX, window_size=window_size,
                                      gap_window_threshold=0.5, gap_column_threshold=0.5, criteria='local', levels=7)

    split_len = int(msa.shape[1] / piece_length) \
        if msa.shape[1] % piece_length == 0 \
        else int(msa.shape[1] / piece_length) + 1

    split_alignment = np.array_split(msa, split_len, axis=1)

    for piece in split_alignment:
        gen_dynamics_vector(piece, dynamics_algorithm)
    print('Number of pieces', int(np.ceil(float(msa.shape[1]) / 5000)))

    #gen_dynamics_vector(msa, dynamics_algorithm=dynamics_algorithm)


    pitch_algorithm = PitchAlgorithm(PitchAlgorithm.WORD_DISTANCES, scale_vector=scale.MajorScale(), n_nucleotides=1)
    durations_algorithm = DurationsAlgorithm(DurationsAlgorithm.FREQUENCIES_DYNAMIC, window_size=1000,
                                             window_duration=500,
                                             n_nucleotides=1)

    similarities_df = pd.DataFrame(np.empty(np.floor(float(msa.shape[1]) / 5000).astype(int),
                                            dtype=[('idx', np.uint8), ('jaccard', np.float), ('tempo', np.float)]))
    j_idx = 0

    k = np.random.choice(np.arange(4, 9), 1)[0]
    n_classes = 3

    print 'K =', k

    for p in range(5000, len(msa[0]), 5000):

        piece_length = 5000 if msa.shape[1] - p >= 5000 else msa.shape[1] - p

        subalignment = np.empty((2, msa.shape[0], piece_length), dtype=np.dtype(str))

        subalignment[0] = msa[:, p-piece_length:p]
        subalignment[1] = msa[:, p:p + piece_length]

        print j_idx
        similarities_df['idx'][j_idx] = j_idx
        similarities_df['jaccard'][j_idx] = calc_jaccard_similarities(subalignment, k=k, inter_alignments=True)['jaccard'][0]
        j_idx += 1

        score = stream.Score()

        score_tempo = tempo.MetronomeMark('adagio', 125)
        score.insert(0, score_tempo)

        pieces.append(score)


    print similarities_df
    j_idx = 0

    classes = np.arange(40, 145, (145 - 40) / n_classes, dtype=np.uint8)
    print classes

    similarities_df = similarities_df.sort_values('jaccard')
    jaccard_indices = np.array_split(similarities_df['jaccard'], n_classes)

    class_size = len(jaccard_indices[0])
    print

    j = 0
    for i in range(0, len(classes)):
        similarities_df['tempo'][j:j + class_size] = classes[i]
        print 'L', len(similarities_df['tempo'][j:j + class_size])
        j += class_size

    print similarities_df.sort_values('idx')
    # jaccard_indices['tempo'] = pd.Series(data=np.zeros(len(jaccard_indices['i'])))"""""