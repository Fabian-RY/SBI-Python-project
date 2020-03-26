#!/usr/bin/env python3
# -*- coding: utf-8 -*-


import tkinter as tk
import tkinter.filedialog

def run_promod():
    pass

def select_directory(entry):
    folder = tk.filedialog.askdirectory(parent=top, initialdir='.', title='Select folder')
    entry.set(folder)
    return None

def select_file(entry):
    folder = tk.filedialog.askopenfilename(parent=top, initialdir='.', title='Select folder')
    entry.set(folder)
    return None

top = tk.Tk()
top.title('Welcome to Promod FC')
#top.geometry('800x600')

## Row 1: Directory with pdbs
name = tk.Label(top, text='Promod FC 0.1')
folder_label = tk.StringVar()
log = tk.Entry(top, textvariable=folder_label, state='disabled')
run = tk.Button(top, text='Select directory', command=lambda:select_directory(folder_label))

## Row 2: Fasta file
fasta_label = tk.StringVar()
fasta_entry = tk.Entry(top, textvariable=fasta_label, state='disabled')
select_fasta_button = tk.Button(top, text='Select directory', command=lambda:select_file(fasta_label))

## Row 3: Threshold

name.grid()
run.grid(column=2, row=0)
log.grid(column=1, row=0)

fasta_entry.grid(column=1, row=1)
select_fasta_button.grid(column=2, row=1)
top.mainloop()
