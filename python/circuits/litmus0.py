import sys, getopt, shutil, os, pprint, networkx as nx

sys.path.append("../eqn")
from eqn_syn import Circuit

sys.path.append("../devices")
from base_devices import ElementType, TwoPortElement
from base_devices import *

from scipy import optimize
from scipy.integrate import ode

import numpy as np

import itertools

import matplotlib.pyplot as plt

import math

import control as ct

from sympy import *

import re as regex

nmos_params = {"KP"     : 120e-6,
               "vth"    : 0.8,
               "lambda" : 0.01,
               "L"      : 2,
               "W"      : 10}

pmos_params = {"KP"     : 40e-6,
               "vth"    : 0.9,
               "lambda" : 0.0125,
               "L"      : 2,
               "W"      : 30}

vgs_1_ref = 0

def add_nmos(g, d, s, KP, vth, l_lambda, L, W, ckt):
    ids = vccs_l1_mosfet()
    ids.set_type(ElementType.current_src)
    ids.set_instance("IDS")
    ids.set_params(KP=KP, \
                   vth=vth, \
                   l_lambda=l_lambda, \
                   L=L, \
                   W=W)

    vgs = current_src()
    vgs.set_type(ElementType.current_src)
    vgs.set_instance("VGS")

    ids_ref = ckt.add_edge(d, s, ids)
    vgs_ref = ckt.add_edge(g, s, vgs)

    ids.set_vgs_ref(vgs_ref=vgs_ref)

def add_pmos(g, d, s, KP, vth, l_lambda, L, W, ckt):
    isd = vccs_l1_mosfet()
    isd.set_type(ElementType.current_src)
    isd.set_instance("IDS")
    isd.set_params(KP=KP, \
                   vth=vth, \
                   l_lambda=l_lambda, \
                   L=L, \
                   W=W)

    vsg = current_src()
    vsg.set_type(ElementType.current_src)
    vsg.set_instance("VGS")

    ids_ref = ckt.add_edge(s, d, isd)
    vgs_ref = ckt.add_edge(s, g, vsg)

    isd.set_vgs_ref(vgs_ref=vgs_ref)

def add_circuit(ckt, rser):

    Vs_val = 1.0
    C2_val = 1e-9
    L1_val = 1e-3
    C1_val = 2e-9
    C3_val = 3e-9
    R2_val = 1e3
    R3_val = 1e-3

    if (rser):
        v_node = 4
    else:
        v_node = 1

    #Add Vs
    Vs = voltage_src(ramp=False, ramp_ddt=1e6)
    Vs.set_type(ElementType.voltage_src)
    Vs.set_value(Vs_val)
    Vs.set_instance("Vs")
    Vs.set_is_input()
    ckt.add_edge(v_node, 0, Vs)

    # Add C2
    C2 = capacitor()
    C2.set_type(ElementType.capacitor)
    C2.set_value(C2_val)
    C2.set_instance("C2")
    ckt.add_edge(1, 3, C2)

    # Add L1
    L1 = inductor()
    L1.set_type(ElementType.inductor)
    L1.set_value(L1_val)
    L1.set_instance("L1")
    ckt.add_edge(1, 2, L1)

    # Add C1
    C1 = capacitor()
    C1.set_type(ElementType.capacitor)
    C1.set_value(C1_val)
    C1.set_instance("C1")
    ckt.add_edge(2, 0, C1)

    # Add C3
    C3 = capacitor()
    C3.set_type(ElementType.capacitor)
    C3.set_value(C3_val)
    C3.set_instance("C3")
    ckt.add_edge(2, 3, C3)

    # Add R2
    R2 = resistor()
    R2.set_type(ElementType.resistor)
    R2.set_value(R2_val)
    R2.set_instance("R2")
    ckt.add_edge(3, 0, R2)

    if (rser):
        R3 = resistor()
        R3.set_type(ElementType.resistor)
        R3.set_value(R3_val)
        R3.set_instance("R3")
        ckt.add_edge(v_node, 1, R3)

def circuit_eqn(x, sys, ckt, vs, t):

    # set Vs
    #vs_ref = ckt.get_ref_from_instance("Vs")
    #ckt.set_value(ref=vs_ref, value=vs)

    # current Im
    I = ckt.get_im(x=x, sys=sys, t=t)

    # voltage Vm
    V = ckt.get_vm(x=x, sys=sys, t=t)

    I = ckt.get_degen()

    # Copy matrices
    qf_num = ckt.qf.copy()
    bf_num = ckt.bf.copy()

    eqn_qf = np.dot(qf_num, I)
    eqn_bf = np.dot(bf_num, V)

    return eqn_qf + eqn_bf

def dypc_litmus0(t, sys, ckt):

    root_start = np.ones(ckt.num_edges)
    vs = 1

    root = optimize.root(circuit_eqn, root_start, args=(sys, ckt, vs, t), tol=1e-7)
    root_start = root.x
    if not root.success:
        print (root.x)
        sys.exit(0)

    return ckt.get_dy(x=root.x)

def ode_solve(ckt):

    num_sys_vars    = ckt.get_num_sys_vars()

    tstep = 10000

    t     = np.linspace(0, 50e-6, tstep)
    x0    = np.zeros (num_sys_vars)

    r = ode(dypc_litmus0).set_integrator('lsoda', method='bdf', atol=1e-9, rtol=1e-9)
    r.set_initial_value(x0, 0.0)
    r.set_f_params(ckt)

    y    = np.empty((tstep, num_sys_vars))
    y[0] = x0

    k = 1
    while r.successful() and k < tstep:
        r.integrate(t[k])

        y[k] = r.y

        k += 1

    return y

def extract_numeric_values(expression):
    numeric_values = regex.findall(r'[-+]?\d+\.\d+', expression)
    numeric_values = [float(value) for value in numeric_values]
    return numeric_values

def separate_numeric_non_numeric(input_string):
    non_numeric_list = regex.findall(r'\b(?<!\d)(?<![a-zA-Z_])[a-zA-Z_]+[a-zA-Z0-9_]*\b', input_string)
    return non_numeric_list

#def extract_numeric_values(expression):
#    numeric_values = regex.findall(r'-?\d+\.\d+|-?\d+(?![\w.])', expression)
#    numeric_values = [float(value) for value in numeric_values]
#    return numeric_values

def map_coefficients_to_variables(expression, numeric_values):
    variables = ['i_L1', 'v_C1', 'v_C2', 'v_C3', 'v_C4']
    coefficient_map = {var: 0.0 for var in variables}

    # Assign the coefficients to respective variables based on index
    for i, var in enumerate(variables):
        coefficient_map[var] = numeric_values[i]

    return coefficient_map

def create_coefficient_vector(expression, numeric_values):
    variables = ['i_L1', 'v_C1', 'v_C2', 'v_C3', 'v_C4']
    coefficient_vector = []

    # Assign the coefficients to a list based on index order
    for i, var in enumerate(variables):
        coefficient_vector.append(numeric_values[i])

    return coefficient_vector

def sym_solve():

    var('i_VS i_C2 v_L1 i_C1 i_C3 i_R2 i_R3 v_C2 v_C1 v_C3 v_VS i_L1 A B C Vaux x u')

    variable_names = ['i_VS', 'i_C2', 'v_L1', 'i_C1', 'i_C3', 'i_R2', 'i_R3']

    # Create a string representation of the list without quotes
    formatted_list = '[' + ', '.join(variable_names) + ']'

    print (formatted_list)

    #Vaux = Matrix([i_VS, i_C2, v_L1, i_C1, i_C3, i_R2, i_R3])
    Vaux = Matrix(variable_names)
    #Vaux = Matrix([formatted_list])
    x    = Matrix([v_C2, v_C1, v_C3, i_L1])

    #print (x)

    u    = Matrix([v_VS])

    #                  +  *  +  +
    A    = Matrix([[1, 0, 0, 0, 0, 0, -1],
                   [0, 1, 0, 0, 0, 0, -1],
                   [0, 0, 1, 0, 0, 0, 0],
                   [0, 0, 0, 1, 0, 1, -1],
                   [0, 0, 0, 0, 1,-1, 1],
                   [0, 0, 0, 0, 0, 1e3,0],
                   [0, 0, 0, 0, 0, 0, 1e-3]])

    #B     = Matrix([[0, 0, 0, 0],
    #                [0, 0, 0, 1],
    #                [-1,0, 1, 0],
    #                [0, 0, 0, 0],
    #                [0, 0, 0,-1],
    #                [0,-1, 1, 0],
    #                [1, 1,-1, 0]])

    B     = Matrix([[ 0,  0,  0,  0],
                    [ 0,  0,  0,  1],
                    [ 0, -1,  1,  0],
                    [ 0,  0,  0,  0],
                    [ 0,  0,  0, -1],
                    [-1,  0,  1,  0],
                    [ 1,  1, -1,  0]] )

    #print (B)

    C     = Matrix([0, 0, 0, 0, 0, 0, -1])

    Ainv  = A.inv()

    i = A.multiply(Vaux) + B.multiply(x) + C.multiply(u)

    j = -A.inv().multiply(B).multiply(x) - A.inv().multiply(C).multiply(u)

    print(type(regex))  # Output the type of the re module
    print(hasattr(regex, 'findall'))  # Check if the findall method is present

    print (j)

    print()
    print()
    print()
    print (j[1])

    string_expression = str(j[1])
    string_expression = regex.sub(r"\s+", "", string_expression)
    numeric_values_list = extract_numeric_values(string_expression)
    coefficients_mapping = create_coefficient_vector(string_expression, numeric_values_list)
    coefficients         = separate_numeric_non_numeric(string_expression)
    print(coefficients_mapping)
    print(coefficients)

def main():

    ckt = Circuit()

    rser = True

    add_circuit(ckt, rser)

    ckt.init_circuit()

    ckt.get_ss()

    sys = ct.ss(ckt.A, ckt.B, ckt.C, ckt.D)
    #mag, phase, omega = ct.freqresp(sys, [0.1, 1., 10.])
    T, yout = ct.step_response(sys, output=0)

    outputs = np.zeros(len(yout[0][0]))

    for i in range(len(outputs)):
        outputs[i] = yout[0][0][i]

    plt.plot(T, outputs)
    plt.show(block=False)

    #exit(0)

    y = ode_solve(ckt)

    fig, ax1 = plt.subplots()
    ax2 = ax1.twinx()

    tstep = 10000
    t     = np.linspace(0, 50e-6, tstep)
    #ax1.plot(T, outputs)
    for sys_var in range(len(y[0,:])):

        inst = str(ckt.get_edge_info_from_sys_var_ref(sys_var).instance)

        if inst[0] == "C":
            l = "$V_{" + inst + "}$"
            ax1.plot(t, y[:,sys_var], label=l)
        else:
            l = "$i_{" + inst + "}$"
            ax2.plot(t, y[:,sys_var], label=l, linestyle="dotted")

    if rser == False:

        VC3   = np.zeros(len(t))

        for idx in range(len(t)):
            for sys_var in range(len(y[0,:])):
                ckt.scb.v[ckt.get_edge_info_from_sys_var_ref(sys_var).ref] = y[idx, sys_var]
            ckt.scb.v[0] = 1.0

            C3 = ckt.get_edge_info(ckt.get_ref_from_instance("C3"))
            C3.get_voltage(ckt.scb)
            VC3[idx] = ckt.scb.v[ckt.get_edge_info(ckt.get_ref_from_instance("C3")).ref]

        l = "$V_{C3}$"
        ax1.plot(t, VC3, label=l)

    ax1.legend()
    ax2.legend(loc="lower right")

    ax1.set_ylabel("$V (V)$")
    ax1.set_xlabel("$t (S)$")

    ax2.set_ylabel("$i (A)$")
    ax2.set_xlabel("$t (S)$")

    if rser == False:
        plt.savefig("litmus0_rser0.png", bbox_inches = 'tight')
    else:
        plt.savefig("litmus0_rser1.svg", bbox_inches = 'tight')

    plt.show()

if __name__ == "__main__":

    main()
