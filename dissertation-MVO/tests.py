from Bio import AlignIO
from music21 import note, converter, duration
from music21 import stream, instrument
from music21.instrument import StringInstrument, WoodwindInstrument, Percussion, PitchedPercussion, UnpitchedPercussion, \
    BrassInstrument
from music21.midi import realtime

"""SCORE_PATH = '/usr/bin/mscore'

n = note.Note("B-2")
half = duration.Duration('half')
n.duration = half

m = note.Note("D2")
m.duration = duration.Duration('quarter')

simple_stream = stream.Stream()

simple_stream.append(n)
simple_stream.append(m)
simple_stream.repeatAppend(n, 2)


# #simple_stream.show(app='/usr/bin/mscore')
simple_stream.show('text')
simple_stream.id = 'My melody'

print simple_stream.highestTime

another_stream = stream.Stream()

nota1 = note.Note("D#5")
nota1.duration = duration.Duration('quarter')

another_stream.append(nota1)

nota2 = note.Note("C#5")
nota2.duration = duration.Duration('quarter')

another_stream.append(nota2)

nota3 = note.Note("C#5")
nota3.duration = duration.Duration('eighth')

another_stream.append(nota3)

simple_stream.append(another_stream)

for thisNote in simple_stream:

    if isinstance(thisNote, stream.Stream):
        for anotherNote in thisNote:
            print anotherNote.step
            print anotherNote.duration
            print anotherNote.pitch.accidental
            print anotherNote.nameWithOctave
    else:
        print thisNote.step
        print thisNote.duration
        print thisNote.pitch.accidental
        print thisNote.nameWithOctave

#####################3
a = note.Note('A4')
a.duration = duration.Duration('half')
a.id = 'A4'

g = note.Note('G4')
g.duration = duration.Duration('quarter')
g.id = 'G4'

d = note.Note('D5')
d.id = 'D5'
d.duration = duration.Duration('quarter')

p = stream.Part()
p.insert(0, instrument.Clarinet())
p.append(a)

p2 = stream.Part()
p2.insert(0, instrument.BassClarinet())
p2.append(g)

score = stream.Score()
score.insert(0, p)
score.insert(1, p2)

import os

file_path = '/home/espinha/Documents/Tese/dissertacao_68794/dissertation-MVO/file.mid'
if os.path.exists(file_path):
    os.unlink(file_path)

#print score.getElementsByClass('Part')[0]
#score.write('midi', fp=file_path)
#score.show(app=SCORE_PATH)

"""""

print PitchedPercussion.__subclasses__()

import os
import sys
import numpy as np

# arr = np.array([1,2,3,4,5,6,7,7,8,9])

"""""
s = stream.Stream()
s.append(n)


sc = stream.Score()
sc.append(s)
print s.beatDuration
"""