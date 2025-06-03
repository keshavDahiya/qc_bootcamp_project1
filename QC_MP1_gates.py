"""
Filename: QC_MP1_gates.py
Author: Keshav Dahiya
Date: 2025-05-30
Description: 
This script contains basic gates needed to implement constrained optimization 
of integer-valued functions of n boolean variables.
The only gates used in the implementation are MCX gates, 1-qubit phase gates and h-gates.
"""

from qiskit.circuit import QuantumCircuit, QuantumRegister, AncillaRegister

import numpy as np

"""
The following function implements the swap gate using 3 cx gates.
"""
def SWAP():
    reg = QuantumRegister(2)
    out = QuantumCircuit(reg, name="SWAP")
    out.cx(reg[0], reg[1])
    out.cx(reg[1], reg[0])
    out.cx(reg[0], reg[1])
    return out.to_gate()

"""
The following function implements the controlled phase gate using 2 cx and 3 1-qubit phase gates.
Input: theta is phase shift
"""
def CP (theta):
    ctr_reg = QuantumRegister(1)
    tgt_reg = QuantumRegister(1)
    out = QuantumCircuit(ctr_reg, tgt_reg, name="CP("+f"{theta/np.pi:.2f}π"+")")
    out.cx(ctr_reg, tgt_reg)
    out.p(-theta/2, tgt_reg)
    out.cx(ctr_reg, tgt_reg)
    out.p(theta/2, ctr_reg)
    out.p(theta/2, tgt_reg)
    return out.to_gate()

"""
The following function implements the multi-controlled phase gate without using any ancillas.
The function relies on recusion. So MCP1(n, theta) will call MCP1(n-1, theta).
In the case of (n-1) control qubits and 1 target qubit, it uses 2(n-1) mcx gates and (2n-1) phase gates.
Inputs: n is the total qubits (the first n-1 being control qubits) and theta is the phase shift.
"""
def MCP1 (n, theta):
    reg = QuantumRegister(n)
    out = QuantumCircuit(reg, name="MCP1("+f"{theta/np.pi:.2f}π"+")")
    if n == 1:
         out.p(theta, reg[0])
         return out.to_gate()
    out.mcx(reg[:-1], reg[-1])
    out.p(-theta/2, reg[-1])
    out.mcx(reg[:-1], reg[-1])
    out.p(theta/2, reg[-1])
    out.compose(MCP1(n-1, theta/2), qubits=reg[:-1], inplace=True)
    return out.to_gate()

"""
The following function implements the multi-controlled phase gate 
using one ancilla, 3 MCX gates and 3 1-qubit phase gates.
Inputs: n is the total qubits (the first n-1 being control qubits) and theta is the phase shift.
"""
def MCP (n, theta):
    ctr_reg = QuantumRegister(n-1)
    tgt_reg = QuantumRegister(1)
    anc_reg = AncillaRegister(1)
    out = QuantumCircuit(ctr_reg, tgt_reg, anc_reg, name="MCP("+f"{theta/np.pi:.2f}π"+")")
    out.mcx(ctr_reg, anc_reg)
    out.mcx(anc_reg, tgt_reg)
    out.p(-theta/2, tgt_reg)
    out.mcx(anc_reg, tgt_reg)
    out.p(theta/2, tgt_reg)
    out.p(theta/2, anc_reg)
    return out.to_gate()

"""
The following function implements the inverse quantum Fourier transform.
It uses the CP and SWAP gates implemented above, and h-gates.
The total gates used for an n-qubit input vector are:
n(n-1)/2 CP gates, n h-gates and n/2 SWAP gates.
Input: n is the number of input qubits.
"""
def IQFT(n):
    reg = QuantumRegister(n)
    out = QuantumCircuit(reg, name="IQFT")
    for i, q in enumerate(reversed(reg), start=1):
        for j, p in enumerate(reversed(reg[n + 1 - i:]), start=1):
            out.compose(CP(-np.pi/(1 << (i-j))), qubits=[q,p], inplace=True)
        out.h(q)
    for q, p in zip(reg[:n >> 1], reversed(reg[n >> 1:])):
        out.compose(SWAP(), qubits=[q,p], inplace=True)
    return out.to_gate()

"""
The following function implements the gate producting a geometric sequence.
For an m-qubit input, the i-th qubit has its phase shifted by theta*(2^i).
It uses m 1-qubit phase gates.
Input: m is the number of input qubits, angle theta is the base phase sift.
"""
def UG (m, theta):
    reg = QuantumRegister(m)
    out = QuantumCircuit(reg, name="UG("+f"{theta/np.pi:.2f}π"+")")
    for i in range(m):
        out.p((1 << i)*theta, reg[i])
    return out.to_gate()

"""
The following function implements a controlled version of the UG gate above.
For an m-qubit input and 1 control qubit, it uses m CP gates implemented above.
Input: m is the number of input qubits, angle theta is the base phase sift.
The first qubit is the control qubit and last m are target qubits.
"""
def CUG (m, theta):
    ctr_reg = QuantumRegister(1)
    tgt_reg = QuantumRegister(m)
    out = QuantumCircuit(ctr_reg, tgt_reg, name="CUG("+f"{theta/np.pi:.2f}π"+")")
    for i in range(m):
        out.compose(CP((1 << i)*theta), qubits=[ctr_reg[0], tgt_reg[i]], inplace=True)
    return out.to_gate()

"""
The following function implements a multi-controlled version of the UG gate above.
For an m-qubit input, it uses m MCP gates implemented above using ancillas.
Input: n is the number of control qubit, m is the number of target qubits, theta is base phase sift.
The first n qubits are the control qubits and last m qubits are the target qubits.
"""
def MCUG (n, m, theta):
    ctr_reg = QuantumRegister(n)
    tgt_reg = QuantumRegister(m)
    anc_reg = AncillaRegister(m)
    out = QuantumCircuit(ctr_reg, tgt_reg, anc_reg, name="MCUG("+f"{theta/np.pi:.2f}π"+")")
    for i in range(m):
        out.compose(MCP(n+1,(1 << i)*theta), qubits=ctr_reg[:]+[tgt_reg[i], anc_reg[i]], inplace=True)
    return out.to_gate()