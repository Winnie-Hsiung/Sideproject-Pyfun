# pip install midiutil

import random
import numpy as np
from midiutil.MidiFile import MIDIFile

class MusicCompositionMachine:
    """A simple music composition machine that generates melodies."""
    
    # Common scales
    SCALES = {
        'C_major': [60, 62, 64, 65, 67, 69, 71, 72],  # C, D, E, F, G, A, B, C
        'C_minor': [60, 62, 63, 65, 67, 68, 70, 72],  # C, D, Eb, F, G, Ab, Bb, C
        'G_major': [55, 57, 59, 60, 62, 64, 66, 67],  # G, A, B, C, D, E, F#, G
        'A_minor': [57, 59, 60, 62, 64, 65, 67, 69],  # A, B, C, D, E, F, G, A
    }
    
    # Common chord progressions (as scale degrees)
    CHORD_PROGRESSIONS = {
        'pop': [1, 5, 6, 4],  # I-V-vi-IV
        'blues': [1, 4, 1, 5, 4, 1],  # I-IV-I-V-IV-I
        'jazz': [2, 5, 1, 6],  # ii-V-I-vi
        'rock': [1, 4, 5, 4]   # I-IV-V-IV
    }
    
    def __init__(self, scale='C_major', tempo=120):
        """Initialize the composition machine with a scale and tempo."""
        self.scale = self.SCALES[scale]
        self.tempo = tempo
        self.midi = MIDIFile(1)  # One track
        self.midi.addTempo(0, 0, tempo)
    
    def generate_melody(self, num_bars=4, beats_per_bar=4, style='random'):
        """Generate a melody with the specified number of bars."""
        notes = []
        durations = []
        
        # Create a rhythm pattern (note durations)
        for _ in range(num_bars):
            if style == 'random':
                # Random rhythmic pattern
                bar_durations = self._generate_random_rhythm(beats_per_bar)
            elif style == 'structured':
                # More structured rhythm
                bar_durations = self._generate_structured_rhythm(beats_per_bar)
            else:
                # Default to simple quarter notes
                bar_durations = [1] * beats_per_bar
            
            durations.extend(bar_durations)
        
        # Create a melodic pattern
        for i in range(len(durations)):
            if style == 'random':
                # Random note from the scale
                note = random.choice(self.scale)
            elif style == 'structured':
                # More structured melodic movement
                if i == 0 or len(notes) == 0:
                    note = self.scale[0]  # Start on the root
                else:
                    # Prefer stepwise motion with occasional jumps
                    prev_index = self.scale.index(notes[-1]) if notes[-1] in self.scale else 0
                    movement = random.choices([-2, -1, 0, 1, 2], weights=[1, 3, 2, 3, 1])[0]
                    new_index = max(0, min(len(self.scale) - 1, prev_index + movement))
                    note = self.scale[new_index]
            else:
                # Simple ascending scale
                note = self.scale[i % len(self.scale)]
            
            notes.append(note)
        
        return notes, durations
    
    def _generate_random_rhythm(self, beats_per_bar):
        """Generate a random rhythm that fits within the beats per bar."""
        durations = []
        remaining_beats = beats_per_bar
        
        while remaining_beats > 0:
            # Possible durations: 0.25 (16th), 0.5 (8th), 1 (quarter), 2 (half)
            possible_durations = [0.25, 0.5, 1, 2]
            # Filter durations that would fit in remaining beats
            valid_durations = [d for d in possible_durations if d <= remaining_beats]
            
            if not valid_durations:
                duration = 0.25  # Default to smallest duration if somehow we get stuck
            else:
                duration = random.choice(valid_durations)
            
            durations.append(duration)
            remaining_beats -= duration
        
        return durations
    
    def _generate_structured_rhythm(self, beats_per_bar):
        """Generate a more musical rhythm pattern."""
        # Common rhythm patterns for different time signatures
        if beats_per_bar == 4:
            patterns = [
                [1, 1, 1, 1],  # All quarter notes
                [2, 2],        # Two half notes
                [1, 0.5, 0.5, 1, 1],  # Quarter, two eighths, two quarters
                [1, 1, 2],     # Two quarters, one half
                [0.5, 0.5, 0.5, 0.5, 2]  # Four eighths, one half
            ]
            return random.choice(patterns)
        else:
            # For other time signatures, default to simple pattern
            return [1] * beats_per_bar
    
    def apply_chord_progression(self, notes, durations, progression_type='pop'):
        """Modify the melody to follow a chord progression."""
        progression = self.CHORD_PROGRESSIONS.get(progression_type, [1, 4, 5, 1])
        
        # Repeat progression as needed to match the melody length
        extended_progression = []
        for i in range(0, len(notes)):
            chord_index = i // 4 % len(progression)  # Change chord every 4 beats
            extended_progression.append(progression[chord_index])
        
        # Adjust notes to fit chord tones when appropriate
        for i in range(len(notes)):
            # Every first beat of a chord, try to use chord tones
            if i % 4 == 0:
                chord_root = extended_progression[i]
                chord_tone = self.scale[(chord_root - 1) % len(self.scale)]
                if random.random() < 0.7:  # 70% chance to use chord tone
                    notes[i] = chord_tone
        
        return notes, durations
    
    def save_midi(self, notes, durations, filename="composition.mid"):
        """Save the notes and durations as a MIDI file."""
        track = 0
        channel = 0
        time = 0  # Start at the beginning
        
        # Add notes to the MIDI file
        for i, (note, duration) in enumerate(zip(notes, durations)):
            self.midi.addNote(track, channel, note, time, duration, 100)  # velocity=100
            time += duration
        
        # Write the MIDI file
        with open(filename, "wb") as output_file:
            self.midi.writeFile(output_file)
        
        return filename
    
    def generate_composition(self, num_bars=8, style='structured', progression='pop', filename="composition.mid"):
        """Generate a complete composition and save it as a MIDI file."""
        notes, durations = self.generate_melody(num_bars, 4, style)
        notes, durations = self.apply_chord_progression(notes, durations, progression)
        return self.save_midi(notes, durations, filename)


# Example usage
if __name__ == "__main__":
    # Create a composition machine with C major scale at 120 BPM
    composer = MusicCompositionMachine(scale='C_major', tempo=120)
    
    # Generate different styles of compositions
    composer.generate_composition(num_bars=8, style='random', progression='pop', filename="random_pop.mid")
    composer.generate_composition(num_bars=8, style='structured', progression='blues', filename="structured_blues.mid")
    composer.generate_composition(num_bars=16, style='structured', progression='jazz', filename="structured_jazz.mid")
    
    print("Music compositions have been created!")
