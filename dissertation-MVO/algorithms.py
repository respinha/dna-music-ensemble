import numpy as np
from music21.instrument import PitchedPercussion
from music21 import scale


# abstract algorithms class
# cannot be instantiated by itself
class Algorithm(dict):

    params = []
    valid_params = []
    valid_algorithms = []

    def __init__(self, algorithm, **kwargs):
        super(Algorithm, self).__init__(**kwargs)

        assert algorithm in self.valid_algorithms
        self['algorithm'] = algorithm

    def __repr__(self):

        ret = self.__str__() + ':\n'

        for p in self.valid_params:
            ret += p
            ret += ','
        ret += '\n'

        ret += 'Items: \n'

        for key, val in self.items():
            ret += key + ': ' + str(val)

        return ret


# pitch algorithms
# supported types: [ dummy, discrete_sync, patternized ]
class DurationsAlgorithm(Algorithm):

    valid_algorithms = ['frequencies_discrete', 'frequencies_dynamic']
    valid_params = ['window_size', 'window_duration', 'n_nucleotides']

    FREQUENCIES_DYNAMIC = 'frequencies_dynamic'
    FREQUENCIES_DISCRETE = 'frequencies_discrete'

    def __init__(self, algorithm='frequencies', **kwargs):
        super(DurationsAlgorithm, self).__init__(algorithm, **kwargs)

        if kwargs:
            for key, value in kwargs.items():
                assert key in self.valid_params, 'Invalid key for ' + self.__str__() + ': ' + key

                if key == 'window_duration':
                    assert isinstance(value, float) or isinstance(value, int)
                else:
                    assert isinstance(value, int)

                self[key] = value

    def __str__(self):
        return 'Durations algorithm'


# pitch algorithms
# supported types: [ word_distances, static_assign ]
class PitchAlgorithm(Algorithm):

    valid_algorithms = ['static_assign', 'word_distances']
    valid_params = ['scale_vector', 'n_nucleotides']

    STATIC_ASSIGN = 'static_assign'
    WORD_DISTANCES = 'word_distances'

    def __init__(self, algorithm, **kwargs):
        super(PitchAlgorithm, self).__init__(algorithm, **kwargs)

        if kwargs:
            for key, value in kwargs.items():
                assert key in self.valid_params, 'Invalid key for ' + self.__str__() + ': ' + key

                if key == 'scale_vector':

                    # assert isinstance(value, list) or isinstance(value, np.ndarray) or is todo: rever estrutura a usar
                    assert isinstance(value, scale.Scale)
                else:
                    assert isinstance(value, int)

                self[key] = value


# dynamics algorithms
# supported types: [ 'gaps' ]
class DynamicsAlgorithm(Algorithm):

    valid_algorithms = ['shannon_index', 'simpson_index']
    valid_params = ['window_size', 'gap_column_threshold', 'gap_window_threshold','criteria', 'levels']

    SHANNON_INDEX = 'shannon_index'
    SIMPSON_INDEX = 'simpsonl_index'

    def __init__(self, algorithm, **kwargs):
        super(DynamicsAlgorithm, self).__init__(algorithm, **kwargs)

        if kwargs:
            for key, value in kwargs.items():
                assert key in self.valid_params, 'Invalid key for ' + self.__str__() + ': ' + key

                if key == 'scale_vector':
                    assert isinstance(value, list) or isinstance(value, np.ndarray)
                elif key == 'criteria':
                    assert isinstance(value, str)
                else:
                    assert isinstance(value, int) or isinstance(value, np.int) or isinstance(value, np.float) or isinstance(value, float)

                self[key] = value

    def __str__(self):
        return 'Dynamics algorithm'


# species clustering
class ClusteringAlgorithm(Algorithm):

    valid_algorithms = ['kmeans', 'hierarchical']
    valid_params = ['n_clusters', 'max_d', 'instrument_pool']

    KMEANS = 'kmeans'
    HIERARCHICAL = 'hierarchical'

    def __init__(self, algorithm, **kwargs):
        super(ClusteringAlgorithm, self).__init__(algorithm)

        if kwargs:
            for key, value in kwargs.items():
                assert key in self.valid_params, 'Invalid key for ' + self.__str__() + ': ' + key

                if key == 'instruments_pool':
                    assert (isinstance(value, np.ndarray) or isinstance(value, list)) and \
                           (all(isinstance(x, PitchedPercussion) for x in value) or all(isinstance(x, str) for x in value))

                elif key == 'n_clusters':
                    assert isinstance(value, int)
                else:
                    assert isinstance(value, float) or isinstance(value, int)

                self[key] = value

    def __str__(self):
        return 'Clustering algorithm'
