# Digital Communications Projects

This repository contains four projects related to digital communications, covering topics such as line coding, matched filters, single-carrier systems, and OFDM simulation.
Each project is implemented in MATLAB and includes analysis, simulation, and visualization of key communication concepts.

---

## 📁 Project 1: Line Coding and Random Process Analysis

### 📌 Objective:
Analyze various line coding schemes and study the statistical properties of the transmitted signals.

### 📋 Tasks:
- Generate ensembles of 500 waveforms for:
  - Unipolar NRZ
  - Polar NRZ
  - Return-to-Zero (RZ)
- Compute:
  - Statistical mean
  - Stationarity
  - Ensemble autocorrelation
  - Time mean and autocorrelation
  - Ergodic property
  - Bandwidth of the transmitted signal

---

## 📁 Project 2: Matched Filters, Correlators, and ISI

### 📌 Objective:
Simulate matched filters and correlators in both noise-free and noisy environments, and analyze inter-symbol interference (ISI) using raised cosine filters.

### 📋 Tasks:
- Implement matched filters and correlators
- Analyze BER vs. Eb/N0 for:
  - Matched filter
  - Non-matched filter
- Compare with theoretical BER
- Study ISI using eye diagrams for different roll-off factors and filter delays

---

## 📁 Project 3: Single-Carrier Communication Systems

### 📌 Objective:
Simulate and analyze the performance of various modulation schemes in an AWGN channel.

### 📋 Tasks:
- Implement baseband communication system with:
  - BPSK, QPSK, 8PSK, BFSK, 16QAM
- Plot BER vs. Eb/N0 for each scheme
- Compare with theoretical BER curves
- Analyze BFSK:
  - Basis functions
  - Baseband equivalent signals
  - Power Spectral Density (PSD)
    
---

## 📁 Project 4: Advanced Communication Systems (OFDM)

### 📌 Objective:
Explore advanced topics including DFT/FFT performance, BER over fading channels, and OFDM system simulation.

### 📋 Tasks:

#### 1. DFT vs. FFT Execution Time:
- Implement DFT
- Compare execution time with FFT

#### 2. BER over Rayleigh Fading Channel:
- Simulate BPSK and 16-QAM
- Compare with and without repetition coding

#### 3. OFDM System Simulation:
- Implement OFDM with coding, interleaving, and cyclic prefix
- Compare BER for:
  - Rayleigh flat fading
  - Frequency-selective fading
- Use BPSK and 16-QAM modulations

---
