"""
Filename: QC_MP1_marker.py
Author: Keshav Dahiya
Date: 2025-05-30
Description: 
This script contains the marker oracle and function encoders needed to implement 
constrained optimization of integer-valued functions of n boolean variables.
"""

from qiskit.circuit import QuantumCircuit, QuantumRegister, AncillaRegister

import numpy as np

from QC_MP1_gates import *

"""
The following function prepares or loads a quantum register with a given number.
Input: num is the given number and in_qubits are the number of qubits to use.
"""
def prep_input (num, in_qubits):
    reg = QuantumRegister(in_qubits)
    out = QuantumCircuit(reg, name="prep_input")
    out.h(reg)
    out.compose(UG(in_qubits, num*np.pi/(1 << (in_qubits-1))), qubits=reversed(reg), inplace=True)
    out.compose(IQFT(in_qubits), qubits=reversed(reg), inplace=True)
    return out.to_gate()

"""
The following function encdes a given function f from F_2^n to Z in a quantum circuit.
Inputs: 1. function is the given function as an array of integer values
        (for example, if f(0)=1, f(1)=2, f(2)=0, f(3)=-3 then function = [1, 2, 0, -3])
        2. out_qubits gives the range of integers for the output of the given function
        (for example, if range is from -8 to 7 then out_qubits=4)
        The output is encoded in 2's complement binary notation.
        3. f_name is the name of the function. For example z=f(x) or z=C(x) etc.
Output: A quantum gate encoding the given function
"""
def function_encoder_anc (function, out_qubits, f_name):
    in_qubits = np.log2(len(function)).astype(int)
    # The following register stores the input x given to a function z=f(x)
    reg_inp = QuantumRegister(in_qubits, name="x")
    # The following register stores the output z of the given function z=f(x)
    reg_out = QuantumRegister(out_qubits, name="z")
    # The following ancilla register is only needed for implementing MCP gate using mcx and 1-qubit phase gates
    reg_anc = AncillaRegister((len(function)-in_qubits-1)*out_qubits, name="a")
    circ = QuantumCircuit(reg_inp, reg_out, reg_anc, name=f_name)
    # The h gates prepare the outut register in a superposition of all states
    circ.h(qubit=reg_out[:])
    theta = np.pi/(1 << out_qubits-1)
    anc = 0 # needed to decide which ancillas to use for MCP gates in an MCUG gate
    for i in range(len(function)):
        bi = format(i, f'0{in_qubits}b')[::-1] #convert index i into a bitstring (can be done using prep_input function as well)
        qbs = [reg_inp[j] for j, x in enumerate(bi) if x == '1'] #load qubits of input register corresponding to a 1 in the bitstring bi
        if  len(qbs) == 0: # no control qubits needed in this case
            circ.compose(UG(out_qubits, function[i]*theta), qubits=reg_out[:], inplace= True)
        elif len(qbs) == 1: # only 1 control qubit needed in this case, so can use CP which does not need any ancilla
            circ.compose(CUG(out_qubits, function[i]*theta), qubits=qbs+reg_out[:], inplace=True)
        else: # otherwise need more than one control qubit and hence MCP gates implemented using ancillas
            circ.compose(MCUG(len(qbs), out_qubits, function[i]*theta), qubits=qbs+reg_out[:]+reg_anc[anc*out_qubits:(anc+1)*out_qubits], inplace=True)
            anc = anc + 1 # update anc to keep track of which ancillas to use for MCP gates
    circ.compose(IQFT(out_qubits), qubits=reg_out[:], inplace=True)
    return circ.to_gate()

"""
The following is the main marker oracle function. _anc indicates that it used ancillas.
The number of anicllas used depend of out_qubits (number of qubits used to encode output).
Precisely, the number of anicllas used are (2^n-n-1)*m for an n-qubit to m-qubit function.
For m-qubit output, the possible range of function encoded is -2^(m-1) to 2^(m-1)-1.
Inputs: 1. f is the objective function (function to be optimized)
        2. t is the threshold value to compare the output of f
        3. C is the contrained function
        4. out_qubits is the number of qubits used to encode the outputes of f and C functions.
Output: A quantum gate which is the marker oracle U. 
        U flips the qubit y[2] in the reg_check[2] if f(x)>t and (C(x)>0 or C(x)=0)
"""
def marker_anc (f, t, C, out_qubits):
    in_qubits = np.log2(len(f)).astype(int)
    reg_inp = QuantumRegister(in_qubits, name="x") #stores input x
    reg_out = QuantumRegister(out_qubits, name="z") #stores function output z
    reg_check = QuantumRegister(3, name="y") #stores the qubits which check the given condition
    #reg_check[0] is flipped if f(x)>t, reg_check[1] is flipped if C(x)\ge 0
    #reg_check[2] is flipped if f(x)>t and C(x)\ge 0
    reg_anc = AncillaRegister((len(f)-in_qubits-1)*out_qubits, name="a") #needed only for MCP implementation
    oracle = QuantumCircuit(reg_inp, reg_out, reg_anc, reg_check, name="marker")
    f[0] = f[0] - t #subtract threshold from f and then check if result is > 0
    oracle.compose(function_encoder_anc(f, out_qubits, "z=f(x)"), qubits=reg_inp[:]+reg_out[:]+reg_anc[:], inplace=True)
    #output is in 2's complement notation so if positive then last qubit would be 0
    oracle.x(reg_out[-1])
    oracle.cx(reg_out[-1], reg_check[0])
    #now ensure that output is strictly greater than 0 by ruling out the possibility of 0 output
    oracle.x(reg_out[:-1])
    oracle.mcx(reg_out[:], reg_check[0])
    #flip to get inp back
    oracle.x(reg_out)
    #take inverse of f encoder
    oracle.compose(function_encoder_anc(f, out_qubits, "z=f(x)").inverse(), qubits=reg_inp[:]+reg_out[:]+reg_anc[:], inplace=True)
    #now check constraint
    oracle.compose(function_encoder_anc(C, out_qubits, "z=C(x)"), qubits=reg_inp[:]+reg_out[:]+reg_anc[:], inplace=True)
    oracle.x(reg_out[-1])
    oracle.cx(reg_out[-1], reg_check[1])
    oracle.x(reg_out[-1])
    #take inverse of C encoder
    oracle.compose(function_encoder_anc(C, out_qubits, "z=C(x)").inverse(), qubits=reg_inp[:]+reg_out[:]+reg_anc[:], inplace=True)
    #if both conditions are satisfied flip the last qubit of y
    oracle.ccx(reg_check[0], reg_check[1], reg_check[2])
    return oracle.to_gate()