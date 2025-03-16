# This is a code for a music creation. It is made for fun. 

#GitHub: https://github.com/nwhitehead/pyfluidsynth

import os
from music21 import stream, note, chord, midi
import fluidsynth
from pydub import AudioSegment

# --- 1. Create a music score with melody and chords, then save as MIDI ---
def create_midi_with_melody_and_chords():
    # Create a music stream
    s = stream.Stream()

    # --- Add melody ---
    melody_notes = [
        note.Note("C4", quarterLength=1),
        note.Note("D4", quarterLength=1),
        note.Note("E4", quarterLength=1),
        note.Note("F4", quarterLength=1),
        note.Note("G4", quarterLength=1),
        note.Note("A4", quarterLength=1),
        note.Note("B4", quarterLength=1),
        note.Note("C5", quarterLength=1)
    ]
    
    # Add melody to the stream
    for melody_note in melody_notes:
        s.append(melody_note)
    
    # --- Add chords ---
    # C Major chord (C4, E4, G4)
    c_major_chord = chord.Chord(["C4", "E4", "G4"])
    c_major_chord.quarterLength = 2
    s.append(c_major_chord)

    # G Major chord (G4, B4, D5)
    g_major_chord = chord.Chord(["G4", "B4", "D5"])
    g_major_chord.quarterLength = 2
    s.append(g_major_chord)

    # F Major chord (F4, A4, C5)
    f_major_chord = chord.Chord(["F4", "A4", "C5"])
    f_major_chord.quarterLength = 2
    s.append(f_major_chord)

    # Save the music score as a MIDI file
    midi_file_path = "output_with_melody_and_chords.mid"
    s.write("midi", fp=midi_file_path)

    print(f"MIDI file with melody and chords saved to {midi_file_path}")
    return midi_file_path

# --- 2. Convert MIDI to audio using FluidSynth ---
def midi_to_audio(midi_file_path):
    soundfont_path = 'path_to_your_soundfont.sf2'  # Set the path to your SoundFont file

    # Initialize FluidSynth
    fs = fluidsynth.Synth()
    fs.sfload(soundfont_path)
    fs.start()

    # Set the output audio file name
    audio_file_path = "output_with_melody_and_chords.wav"

    # Convert MIDI to audio (WAV)
    fs.midi_to_audio(midi_file_path, audio_file_path)

    print(f"Audio with melody and chords saved to {audio_file_path}")
    return audio_file_path

# --- 3. Convert audio format using Pydub ---
def convert_audio_format(input_audio_path, output_audio_path):
    # Use Pydub to load the audio file
    audio = AudioSegment.from_wav(input_audio_path)

    # Convert the audio to MP3 format
    audio.export(output_audio_path, format="mp3")
    print(f"Audio converted to MP3: {output_audio_path}")

# --- 4. Complete process execution ---
def main():
    # 1. Create the music score with melody and chords, then save as a MIDI file
    midi_file = create_midi_with_melody_and_chords()

    # 2. Use FluidSynth to convert the MIDI to audio
    audio_file = midi_to_audio(midi_file)

    # 3. Optional: Convert the WAV audio to MP3 format
    convert_audio_format(audio_file, "output_with_melody_and_chords.mp3")

    # Clean up the intermediate files
    os.remove(midi_file)  # Delete the generated MIDI file
    os.remove(audio_file)  # Delete the generated WAV file

    print("Process complete!")

# Run the main process
if __name__ == "__main__":
    main()
