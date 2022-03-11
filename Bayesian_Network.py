from numpy.core.function_base import linspace
from pgmpy.inference.ExactInference import BeliefPropagation
from pgmpy.inference.mplp import Mplp
from pgmpy.models import BayesianModel
from pgmpy.factors.discrete.CPD import TabularCPD
from pgmpy.inference import VariableElimination
from pgmpy.inference.CausalInference import CausalInference
import time

import pandas as pd
from pgmpy.estimators import BayesianEstimator
import matplotlib.pyplot as plt
from pylab import *
import numpy as np

'''
Node:
    ROOT>>>
    
    short_rad           线径小       
    high_load           负载高       
    humidity            潮湿         
    corrosion           腐蚀         
    cut                 割损         
    biten_by_mouse      老鼠咬       
    term_corrosion      接线端子锈蚀  
    term_loosen         接触端子松动  
    circuit             线圈         
    work_vol            工作电压     
    light_surf_temp     灯具镇流器外壳温度
    elec_heat_tool      电热器具

    INTERMEDIATE>>>

    elec_fire           电器火灾
    circuit_road        线路
    temp                温度
    temp_rise           温升
    con_temp            连接点温度
    side_temp           分支节点温度
    line_term_rise      接线端子温升值
    bad_connection      接触不良
    circuit_temp_rise   线路温升值
    overload            过负荷
    bad_env             环境不良
    broken              破损
    current             电流
    circuit_insu_red    线路绝缘降低
    short_circuit       短路
    elec_leakage        漏电
    elec_spark          电火花
    rect_equip          整流设备
    temp_rise_1         温升1
    low_vol_equip       低压用电设备
'''
def main():
    BN = BayesianModel()


    #添加节点
    BN.add_nodes_from(nodes=['elec_fire', 'circuit_road', 'temp', 'temp_rise', 'con_temp', 'side_temp', 'line_term_rise'
    , 'bad_connection', 'term_corrosion', 'term_loosen', 'circuit_temp_rise', 'overload', 'bad_env', 'broken',
    'current', 'short_rad', 'high_load', 'humidity', 'corrosion', 'cut', 'biten_by_mouse', 'circuit_insu_red',
    'short_circuit', 'elec_leakage', 'elec_spark', 'light_surf_temp', 'elec_heat_tool', 'rect_equip', 'work_vol',
    'temp_rise_1', 'circuit', 'low_vol_equip'])

    #添加节点之间的因果关系
    BN.add_edges_from(ebunch=[('circuit_road','elec_fire'), ('low_vol_equip','elec_fire'), 
    ('temp_rise', 'circuit_road'), ('temp', 'circuit_road'), ('elec_spark', 'circuit_road'), ('light_surf_temp', 'low_vol_equip'),
    ('elec_heat_tool', 'low_vol_equip'), ('rect_equip', 'low_vol_equip'), ('con_temp', 'temp'), ('side_temp', 'temp'),
    ('line_term_rise', 'temp_rise'), ('circuit_temp_rise', 'temp_rise'), ('elec_leakage', 'elec_spark'), ('short_circuit', 'elec_spark'),
    ('work_vol', 'rect_equip'), ('temp_rise_1', 'rect_equip'), ('circuit', 'temp_rise_1'), 
    ('bad_connection', 'con_temp'), ('bad_connection', 'side_temp'), ('bad_connection', 'line_term_rise'), 
    ('overload', 'circuit_temp_rise'), ('circuit_insu_red', 'short_circuit'), ('circuit_insu_red', 'elec_leakage'), 
    ('term_loosen', 'bad_connection'), ('term_corrosion', 'bad_connection'), ('overload', 'circuit_insu_red'),
    ('bad_env', 'circuit_insu_red'), ('broken', 'circuit_insu_red'), ('current', 'overload'), 
    ('humidity', 'bad_env'), ('corrosion', 'bad_env'), ('cut', 'broken'), ('biten_by_mouse', 'broken'), 
    ('short_rad', 'current'), ('high_load', 'current')])

    #定义节点的条件概率分布
    #根节点

    short_rad_cpd = TabularCPD(variable = 'short_rad', variable_card = 2, values=[[0.9],[0.1]])

    high_load_cpd = TabularCPD(variable = 'high_load', variable_card = 2, values=[[0.9],[0.1]])

    humidity_cpd  = TabularCPD(variable = 'humidity', variable_card = 2, values=[[0.9],[0.1]])

    corrosion_cpd = TabularCPD(variable = 'corrosion', variable_card = 2, values=[[0.9],[0.1]])

    cut_cpd  = TabularCPD(variable = 'cut', variable_card = 2, values=[[0.9],[0.1]])

    biten_by_mouse_cpd = TabularCPD(variable = 'biten_by_mouse', variable_card = 2, values=[[0.9],[0.1]])

    term_corrosion_cpd = TabularCPD(variable = 'term_corrosion', variable_card = 2, values=[[0.9],[0.1]])

    term_loosen_cpd = TabularCPD(variable = 'term_loosen', variable_card = 2, values=[[0.9],[0.1]])

    circuit_cpd = TabularCPD(variable = 'circuit', variable_card = 2, values=[[0.9],[0.1]])

    work_vol_cpd = TabularCPD(variable = 'work_vol', variable_card = 2, values=[[0.9],[0.1]])

    light_surf_temp_cpd = TabularCPD(variable = 'light_surf_temp', variable_card = 2, values=[[0.9],[0.1]]) 

    elec_heat_tool_cpd = TabularCPD(variable = 'elec_heat_tool', variable_card = 2, values=[[0.9],[0.1]])    


    #非根节点
    elec_fire_cpd = TabularCPD(variable = 'elec_fire', variable_card = 2, values=[[0.999, 0.71, 0.06, 0.05],
                                                                                [0.001, 0.29, 0.94, 0.95]],
                                                                                    evidence=['circuit_road', 'low_vol_equip'], 
                                                                                    evidence_card=[2, 2])   
    circuit_road_cpd = TabularCPD(variable = 'circuit_road', variable_card = 2, values=[[0.99, 0.80, 0.71, 0.61, 0.3, 0.2, 0.05, 0.01],
                                                                                        [0.01, 0.2, 0.29, 0.39, 0.7, 0.8, 0.95, 0.99]], 
                                                                                        evidence=['elec_spark', 'temp', 'temp_rise'], 
                                                                                        evidence_card=[2, 2, 2])
    print(circuit_road_cpd)       
    temp_cpd = TabularCPD(variable = 'temp', variable_card = 2, values=[[0.99, 0.71, 0.06, 0.05],
                                                                        [0.01, 0.29, 0.94, 0.95]], 
                                                                        evidence=['con_temp', 'side_temp'], 
                                                                        evidence_card=[2, 2])
              
    temp_rise_cpd = TabularCPD(variable = 'temp_rise', variable_card = 2, values=[[0.99, 0.81, 0.06, 0.01],
                                                                                    [0.01, 0.19, 0.94, 0.99]], 
                                                                                    evidence=['line_term_rise', 'circuit_temp_rise'], 
                                                                                    evidence_card=[2, 2])
    #temp_rise_cpd = TabularCPD(variable = 'temp_rise', variable_card = 2, values=[[1, 1, 1, 1],
    #                                                                                [0, 0, 0, 0]], 
    #                                                                                evidence=['line_term_rise', 'circuit_temp_rise'], 
    #                                                                               evidence_card=[2, 2])
    con_temp_cpd = TabularCPD(variable = 'con_temp', variable_card = 2, values=[[0.99, 0.1],
                                                                                [0.01, 0.9]], 
                                                                                evidence=['bad_connection'], 
                                                                                evidence_card=[2])           
    side_temp_cpd = TabularCPD(variable = 'side_temp', variable_card = 2, values=[[0.99, 0.1],
                                                                                    [0.01, 0.9]], 
                                                                                    evidence=['bad_connection'], 
                                                                                    evidence_card=[2])          
    line_term_rise_cpd = TabularCPD(variable = 'line_term_rise', variable_card = 2, values=[[0.99, 0.1],
                                                                                            [0.01, 0.9]], 
                                                                                            evidence=['bad_connection'], 
                                                                                            evidence_card=[2])    
    bad_connection_cpd = TabularCPD(variable = 'bad_connection', variable_card = 2, values=[[0.99, 0.71, 0.06, 0.05],
                                                                                            [0.01, 0.29, 0.94, 0.95]],
                                                                                            evidence=['term_loosen', 'term_corrosion'], 
                                                                                            evidence_card=[2, 2])       
    circuit_temp_rise_cpd = TabularCPD(variable = 'circuit_temp_rise', variable_card = 2, values=[[0.99, 0.01],
                                                                                                [0.01, 0.99]],
                                                                                                evidence=['overload'], 
                                                                                                evidence_card=[2])  
    overload_cpd = TabularCPD(variable = 'overload', variable_card = 2, values=[[0.99, 0.01],
                                                                                [0.01, 0.99]],
                                                                                evidence=['current'], 
                                                                                evidence_card=[2])           
    bad_env_cpd = TabularCPD(variable = 'bad_env', variable_card = 2, values=[[0.99, 0.71, 0.06, 0.05],
                                                                                [0.01, 0.29, 0.94, 0.95]], 
                                                                                evidence=['humidity', 'corrosion'], 
                                                                                evidence_card=[2, 2])           
    broken_cpd = TabularCPD(variable = 'broken', variable_card = 2, values=[[0.99, 0.71, 0.06, 0.05],
                                                                            [0.01, 0.29, 0.94, 0.95]], 
                                                                            evidence=['biten_by_mouse', 'cut'], 
                                                                            evidence_card=[2, 2])             
    current_cpd = TabularCPD(variable = 'current', variable_card = 2, values=[[0.99, 0.81, 0.06, 0.01],
                                                                                [0.01, 0.19, 0.94, 0.99]], 
                                                                                evidence=['short_rad', 'high_load'], 
                                                                                evidence_card=[2, 2])           
    circuit_insu_red_cpd = TabularCPD(variable = 'circuit_insu_red', variable_card = 2, values=[[0.99, 0.80, 0.71, 0.61, 0.3, 0.2, 0.05, 0.01],
                                                                                        [0.01, 0.2, 0.29, 0.39, 0.7, 0.8, 0.95, 0.99]], 
                                                                                        evidence=['overload', 'bad_env', 'broken'], 
                                                                                        evidence_card=[2, 2, 2])  
    short_circuit_cpd = TabularCPD(variable = 'short_circuit', variable_card = 2, values=[[0.99, 0.1],
                                                                                        [0.01, 0.9]], 
                                                                                        evidence=['circuit_insu_red'], 
                                                                                        evidence_card=[2])     
    elec_leakage_cpd = TabularCPD(variable = 'elec_leakage', variable_card = 2, values=[[0.99, 0.1],
                                                                                        [0.01, 0.9]], 
                                                                                        evidence=['circuit_insu_red'], 
                                                                                        evidence_card=[2])     
    elec_spark_cpd = TabularCPD(variable = 'elec_spark', variable_card = 2, values=[[0.99, 0.71, 0.06, 0.05],
                                                                            [0.01, 0.29, 0.94, 0.95]], 
                                                                            evidence=['elec_leakage', 'short_circuit'], 
                                                                            evidence_card=[2, 2])          
    rect_equip_cpd = TabularCPD(variable = 'rect_equip', variable_card = 2, values=[[0.99, 0.71, 0.06, 0.05],
                                                                            [0.01, 0.29, 0.94, 0.95]], 
                                                                                    evidence=['work_vol', 'temp_rise_1'], 
                                                                                    evidence_card=[2, 2])                
    temp_rise_1_cpd = TabularCPD(variable = 'temp_rise_1', variable_card = 2, values=[[0.99, 0.1],
                                                                                    [0.01, 0.9]], 
                                                                                    evidence=['circuit'], 
                                                                                    evidence_card=[2])              
    low_vol_equip_cpd = TabularCPD(variable = 'low_vol_equip', variable_card = 2, values=[[0.99, 0.80, 0.71, 0.61, 0.3, 0.2, 0.05, 0.01],
                                                                                        [0.01, 0.2, 0.29, 0.39, 0.7, 0.8, 0.95, 0.99]], 
                                                                                        evidence=['light_surf_temp', 'elec_heat_tool', 'rect_equip'], 
                                                                                        evidence_card=[2, 2, 2])

    #将条件概率表加入模型
    BN.add_cpds(elec_fire_cpd, circuit_road_cpd, temp_cpd, temp_rise_cpd, con_temp_cpd, side_temp_cpd, 
    line_term_rise_cpd, bad_connection_cpd, term_corrosion_cpd, term_loosen_cpd, circuit_temp_rise_cpd, 
    overload_cpd, bad_env_cpd, broken_cpd, current_cpd, short_rad_cpd, high_load_cpd, humidity_cpd, corrosion_cpd,
    cut_cpd, biten_by_mouse_cpd, circuit_insu_red_cpd, short_circuit_cpd, elec_leakage_cpd, elec_spark_cpd, 
    light_surf_temp_cpd, elec_heat_tool_cpd, rect_equip_cpd, work_vol_cpd, temp_rise_1_cpd, circuit_cpd, 
    low_vol_equip_cpd)

    #检查模型
    print(BN.check_model())
    time_start2 = time.time()
    model_inference2 = BeliefPropagation(BN)
    data = pd.DataFrame(np.random.randint(low = 0, high = 2, size = (1000, size(BN.nodes()))), 
    columns = ['elec_fire', 'circuit_road', 'temp', 'temp_rise', 'con_temp', 'side_temp', 'line_term_rise'
    , 'bad_connection', 'term_corrosion', 'term_loosen', 'circuit_temp_rise', 'overload', 'bad_env', 'broken',
    'current', 'short_rad', 'high_load', 'humidity', 'corrosion', 'cut', 'biten_by_mouse', 'circuit_insu_red',
    'short_circuit', 'elec_leakage', 'elec_spark', 'light_surf_temp', 'elec_heat_tool', 'rect_equip', 'work_vol',
    'temp_rise_1', 'circuit', 'low_vol_equip'])




    #变量消元
    model_inference1 = VariableElimination(BN)
    test = model_inference1.query(variables={'elec_fire'}, 
                                    evidence={'current': 1},
                                    show_progress=False)
    print(test)

    pic(model_inference1)

'''
def pic_plot(factor):
    t = linspace(0, 10, 200)
    
    s = [0 if (i<5 and i>=0) else factor for i in t]
    
    fig, ax = plt.subplots()
    ax.plot(t, s)
    ax.set(xlabel = 'time(s)', ylabel = 'probability', title='Probability of occuring fire')
    ax.grid()
    plt.show()
'''

def pic(inference):
    
    sensor_temperature = 25
    temp_prob = []
    elec_fire_prop = []
    elec_fire_prop_cur = []

    i_cur = 0
    i_temp = 0

    threshold_temp = 100
    time_t = linspace(0, 10, 20)
    temperature_t = linspace(25, 150, 20)
    current_t = linspace (0, 1.5, 20)
    

    temp = [sensor_temperature if(i>=0 and i<=5) else (20*i - 75 + np.random.randint(low = 15, size=1)) for i in time_t] # 原始数据
    sig_array = [sigmoid(j, 100, 1/3) for j in temperature_t]
    current = [0.5 if(k>=0 and k<=5) else (0.1*k + 0.005* np.random.randint(low = 10, size=1)) for k in time_t]
    sig_array_current = [sigmoid(l, 0.8, 20) for l in current_t]

    cur_prob = []
    for cur in current:
        global is_high_current
        i_cur += 1
        sigmoid_index_cur = sigmoid(cur, 0.8, 20)
        cur_prob.append(sigmoid_index_cur)
        if sigmoid_index_cur >= 0.5:
            is_high_current = 1
        else:
            is_high_current = 0
        
        if i_cur >= 10:
            is_high_load_cur = 1
        else:
            is_high_load_cur = 0

        current_prob = inference.query(variables={'elec_fire'}, 
                                    evidence={'current': is_high_current, 'high_load':is_high_load_cur},
                                    show_progress=False)
        factor = current_prob.values[1]
        factor = float(factor)
        elec_fire_prop_cur.append(factor)
    

    for tep in temp:
        if cur_prob[i_temp] >= 0.5: 
            is_high_current_temp = 1 
        else: 
            is_high_current_temp = 0
        i_temp += 1
        sigmoid_index = sigmoid(tep, threshold_temp, 1/3) # Sigmoid转化后温度值
        
        temp_prob.append(sigmoid_index)
        if sigmoid_index >= 0.5:
            is_high_temp = 1
        else:
            is_high_temp = 0
       
        if i_temp >= 10:
            is_high_load = 1
        else:
            is_high_load = 0
        
        
        fire_prop = inference.query(variables={'elec_fire'}, 
                                    evidence={'temp': is_high_temp, 'high_load':is_high_load, 'current': is_high_current_temp},
                                    show_progress=False)

        factor = fire_prop.values[1]
        factor = float(factor)
        elec_fire_prop.append(factor) # 以温度为证据，推理后火灾发生概率

 # 画图

    fig, ax = plt.subplots()
    ax.plot(time_t, temp_prob, label = 'Probability', 
            markersize=6, marker='o',markeredgecolor = 'k', 
            markerfacecolor = 'y')
    ax.set(xlabel = 'Time(min)', 
           ylabel = 'Probability', 
           title = 'Probability reflection of temperature by sigmoid function')
    plt.annotate("Temperature starts to rise", (5,temp_prob[5]), xycoords='data',
             xytext=(2, 0.45), 
             arrowprops=dict(arrowstyle='->'))
    
    plt.legend()
    plt.show()
    
    fig, ax1 = plt.subplots()
    ax1.plot(time_t, elec_fire_prop, c = 'r', label = 'Probability of elec_fire',
            markersize=6, marker='o',markeredgecolor = 'k', 
            markerfacecolor = 'y')
    ax1.set(xlabel = 'Pime(min)', 
            ylabel = 'Probability', 
            title = 'Probability of occuring fire when temp rises')
    plt.legend()
    plt.show()


    fig, ax2 = plt.subplots()
    ax2.plot(time_t, temp, label = 'Temperature', 
            markersize=6, marker='o',markeredgecolor = 'k', 
            markerfacecolor = 'y')
    ax2.set(xlabel = 'Time(min)', 
           ylabel = 'Temperature(°C)', 
           title = 'Temperature Rise Curve')
    
    plt.legend()
    plt.show()



    fig, ax3 = plt.subplots()
    ax3.plot(temperature_t, sig_array, label = 'Sigmoid Function')
    ax3.set(xlabel = 'Temperature(°C)',
            ylabel = 'Probability',
            title = 'Sigmoid function curve of temperature')
    plt.legend()
    plt.show()

    fig, ax4 = plt.subplots()
    ax4.plot(time_t, elec_fire_prop_cur, label = 'Prob',
            markersize=6, marker='o',markeredgecolor = 'k', 
            markerfacecolor = 'y')
    ax4.set(xlabel = 'Time(min)',
            ylabel = 'Probability',
            title = 'Probability of occuring fire when the data is missing')
    plt.legend()
    plt.show()

    fig, ax5 = plt.subplots()
    ax5.plot(time_t, current, label = 'Prob',
            markersize=6, marker='o',markeredgecolor = 'k', 
            markerfacecolor = 'y')
    ax5.set(xlabel = 'Time(min)',
            ylabel = 'Current(A)',
            title = 'Current rising curves')
    plt.legend()
    plt.show()

    ig, ax6 = plt.subplots()
    ax6.plot(current_t, sig_array_current, label = 'Prob')
    ax6.set(xlabel = 'Cuurent(A)',
            ylabel = 'Probability',
            title = 'Sigmoid function curve of current')
    plt.legend()
    plt.show()

    ig, ax7 = plt.subplots()
    ax7.plot(time_t, cur_prob, label = 'Prob')
    ax7.set(xlabel = 'Cuurent(A)',
            ylabel = 'Probability',
            title = 'Current probability')
    plt.legend()
    plt.show()
    


def sigmoid(x, c, k):
    return 1 / (1 + np.exp(k*(-x + c)))
    


if __name__ == "__main__":
    main()