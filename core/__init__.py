from bio import gen_alignment, gen_random_seqs
from music import gen_song

import numpy as np


def load_seq_config(filepath):
    from os import path
    assert path.isfile(filepath), 'Invalid JSON path ' + filepath

    import json
    data = json.load(open(filepath, 'r'))

    assert 'length' in data, 'Length not specified'

    len = data[u'length']
    del data[u'length']

    seqs = np.zeros((np.alen(data.keys()), len), dtype="S1")
    i = 0
    for k in data.keys():

        seq = data[k]

        for elem in seq:
            symbol = elem['symbol']
            idx = elem['idx']

            idx = [int(_idx) for _idx in idx.split(':')]
            assert isinstance(idx, list), 'Invalid data type for boundaries ' + str(type(idx))
            assert np.alen(idx) == 2 and idx[1] > idx[0] and idx[1] <= len, \
                'Invalid boundaries for idx: ' + str(idx[0]) + ', ' + str(idx[1])

            idx = [int(_idx) for _idx in idx]
            if symbol != 'R':
                seqs[i][idx[0]:idx[1]] = symbol
            else:
                seqs[i][idx[0]:idx[1]] = np.random.choice(['a','g','c','t','-'], idx[1] - idx[0])

            print idx[0], idx[1], seqs[i][idx[0]:idx[1]]
        i += 1

    return seqs