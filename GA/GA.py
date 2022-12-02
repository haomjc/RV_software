# -*- coding: utf-8 -*-
from PySide2.QtCore import Slot
from PySide2.QtWidgets import QDialog, QMessageBox, QTableWidgetItem
import matplotlib.pyplot as plt
# 导入matplotlib模块并使用Qt5Agg
'''
import matplotlib
matplotlib.use('Qt5Agg')
# 使用 matplotlib中的FigureCanvas (在使用 Qt5 Backends中 FigureCanvas继承自QtWidgets.QWidget)
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas

import mpl_toolkits.mplot3d as mp3d
'''
from mpl_toolkits.mplot3d import Axes3D  # 空间三维画图
from GA.Ui_GA import Ui_Dialog_GA
import geatpy as ga # 导入geatpy库
import numpy as np
#================= 导入自定义的进化算法模板===========================
from GA.nsga2_usr import nsga2
from GA.awGA_usr import awGA
from GA.rwGA_usr import rwGA
from GA.q_sorted_usr import q_sorted
from GA.q_sorted_new_usr import q_sorted_new

class Dialog_GA(QDialog, Ui_Dialog_GA):

    def __init__(self, parent=None):
        super(Dialog_GA, self).__init__(parent)
        self.setupUi(self)
        self.radioButton_q_sorted = False
        self.radioButton_awGA = False
        self.radioButton_rwGA = False
        self.radioButton_q_sorted_new =False
        self.radioButton_NSGA_II = False
        '''
        # 几个QWidgets
        self.figure = mp3d
        self.canvas = FigureCanvas(self.figure)
        self.verticalLayout.addWidget(self.canvas)
        '''
    @Slot()
    def on_radioButton_q_sorted_clicked(self):
        
        self.radioButton_q_sorted = True
        self.radioButton_awGA = False
        self.radioButton_rwGA = False
        self.radioButton_q_sorted_new =False
        self.radioButton_NSGA_II = False    
        
    @Slot()
    def on_radioButton_awGA_clicked(self):

        self.radioButton_q_sorted = False
        self.radioButton_awGA = True
        self.radioButton_rwGA = False
        self.radioButton_q_sorted_new =False
        self.radioButton_NSGA_II = False
    
    @Slot()
    def on_radioButton_rwGA_clicked(self):

        self.radioButton_q_sorted = False
        self.radioButton_awGA = False
        self.radioButton_rwGA = True
        self.radioButton_q_sorted_new =False
        self.radioButton_NSGA_II = False
    
    @Slot()
    def on_radioButton_q_sorted_new_clicked(self):

        self.radioButton_q_sorted = False
        self.radioButton_awGA = False
        self.radioButton_rwGA = False
        self.radioButton_q_sorted_new =True
        self.radioButton_NSGA_II = False
    
    @Slot()
    def on_radioButton_NSGA_II_clicked(self):

        self.radioButton_q_sorted = False
        self.radioButton_awGA = False
        self.radioButton_rwGA = False
        self.radioButton_q_sorted_new =False
        self.radioButton_NSGA_II = True
    
    @Slot()
    def on_pushButton_design_range_clicked(self):

        self.P = self.doubleSpinBox_P.value()
        self.n = self.doubleSpinBox_n.value()
        self.n_p = self.spinBox_n_p.value()
        self.T_2 = self.doubleSpinBox_T_2.value()    #单位：N/mm    
        
        D_z_min = 2.1*self.T_2**(1/3)        
        D_z_max = 2.6*self.T_2**(1/3)       
        
        self.doubleSpinBox_L_max.setValue(0.075*D_z_max)
        
        self.doubleSpinBox_D_z_min.setValue(D_z_min)
        self.doubleSpinBox_D_z_max.setValue(D_z_max)
        
        self.doubleSpinBox_B_min.setValue(0.075*D_z_min)
        self.doubleSpinBox_B_max.setValue(0.075*D_z_max)
        
    @Slot()
    def on_pushButton_optimal_calc_clicked(self):
        
        # 获取函数接口地址
        AIM_M = __import__('GA.aimfuc')
        PUN_M = __import__('GA.punishing')
        #=============================================输入值=============================================
        
        K_1_min = self.doubleSpinBox_K_1_min.value()
        K_1_max = self.doubleSpinBox_K_1_max.value()
        
        D_r_min = self.doubleSpinBox_D_r_min.value()
        D_r_max = self.doubleSpinBox_D_r_max.value()
        
        z_2_min = self.spinBox_z2_min.value()
        z_2_max = self.spinBox_z2_max.value()
        
        z_g_min = self.spinBox_z_g_min.value()
        z_g_max = self.spinBox_z_g_max.value()
        
        m_min = self.spinBox_m_min.value()
        m_max = self.spinBox_m_max.value()
        
        b_min = self.doubleSpinBox_b_min.value()
        b_max  = self.doubleSpinBox_b_max.value()
        
        z_1_min = self.spinBox_z1_min.value()
        z_1_max = self.spinBox_z1_max.value()        
        
        L_min = self.doubleSpinBox_L_min.value()
        L_max = self.doubleSpinBox_L_max.value()
        
        D_z_min = self.doubleSpinBox_D_z_min.value()
        D_z_max = self.doubleSpinBox_D_z_max.value()
        
        Z_min  = self.spinBox_Z_min.value()
        Z_max  = self.spinBox_Z_max.value()
        
        d_z_min = self.doubleSpinBox_d_z_min.value()
        d_z_max = self.doubleSpinBox_d_z_max.value()
        
        D_m_min = self.doubleSpinBox_D_m_min.value()
        D_m_max = self.doubleSpinBox_D_m_max.value()   
   
        B_min = self.doubleSpinBox_B_min.value()
        B_max = self.doubleSpinBox_B_max.value()
        #============================变量设置============================
        b = [b_min, b_max]     #行星轮模数与齿宽的关系
    
        D_z = [D_z_min, D_z_max]
        
        d_z = [ d_z_min, d_z_max]
        B = [B_min, B_max]  #摆线轮宽度
        K_1 = [K_1_min,K_1_max]
        
        D_m = [D_m_min, D_m_max]
        D_r = [D_r_min, D_r_max]
        
        L = [L_min, L_max]
        
        z_1 = [z_1_min, z_1_max]   #不根切条件
        z_2 = [z_2_min, z_2_max]
        m_4 = [m_min*4, m_max*4]
        z_g_2 = [z_g_min//2, z_g_max//2]
        Z = [Z_min,Z_max]  #滚子数最小值为5
        
        b1 = [1, 1] # 自变量是否包含上下界
        b2 = [1, 1] 
        b3 = [1, 1] 
        b4 = [1, 1] 
        b5 = [1, 1] 
        b6 = [1, 1] 
        b7 = [1, 1] 
        b8 = [1, 1]
        b9 = [1, 1] 
        b10 = [1, 1] 
        b11 = [1, 1] 
        b12 = [1, 1] 
        b13 = [1, 1]
        
        ranges = np.vstack([b, D_z, d_z, B, K_1,D_m, D_r, L, z_1, z_2, m_4, z_g_2, Z]).T  # 生成自变量的范围矩阵
        borders = np.vstack([b1, b2, b3, b4, b5, b6, b7, b8, b9, b10, b11, b12, b13]).T   # 生成自变量的边界矩阵（1表示变量的区间是闭区间）
        precisions = [2, 2, 2, 2, 4, 2, 2, 2, 0, 0, 0, 0, 0]       # 根据crtfld的函数特性，这里需要设置精度为任意正值，否则在生成区域描述器时会默认为整数编码，并对变量范围作出一定调整
        
        #========================遗传算法参数设置=========================
        
        NIND = self.spinBox_NIND.value()           # 种群规模
        MAXGEN = self.spinBox_MAXGEN.value()            # 最大遗传代数
        MAXSIZE = self.spinBox_pareto_max.value()
        GGAP = self.doubleSpinBox_GGAP.value()             # 代沟：子代与父代的重复率为(1-GGAP)  ，为了避免父子两代合并后种群数量爆炸，要让代沟为0.5
        selectStyle = 'tour'     # 遗传算法的选择方式设为"rws"——轮盘赌选择
        recombinStyle = 'xovdp'  # 遗传算法的重组方式，设为两点交叉
        problem = 'M'  #混合型编码
        R_num = 8    #实数变量数量
        
        recopt = self.doubleSpinBox_recopt.value()           # 交叉概率
        pm = self.doubleSpinBox_pm.value()                 # 变异概率
        SUBPOP = 1               # 设置种群数为1
        maxormin = 1             # 设置标记表明这是最小化目标
        
        FieldDR = ga.crtfld(ranges, borders, precisions) # 生成区域描述器
        #=======================调用种群进化模板===================

        if self.radioButton_NSGA_II == True:   #效果好
            [ObjV, NDSet, NDSetObjV, times] = nsga2.moea_nsga2(self, AIM_M, 'aimfuc', PUN_M, 'punishing', FieldDR,problem,R_num,maxormin,MAXGEN, MAXSIZE, NIND,SUBPOP,GGAP,selectStyle, recombinStyle,recopt,pm, distribute = False)
        elif self.radioButton_awGA == True: #效果差
            [ObjV, NDSet, NDSetObjV, times] = awGA.moea_awGA(self,AIM_M, 'aimfuc', PUN_M, 'punishing', FieldDR,problem, R_num,maxormin,MAXGEN, MAXSIZE, NIND,SUBPOP,GGAP,selectStyle, recombinStyle,recopt,pm, distribute = False)
        elif self.radioButton_rwGA == True:   #效果差
            [ObjV, NDSet, NDSetObjV, times] = rwGA.moea_rwGA(self,AIM_M, 'aimfuc', PUN_M, 'punishing', FieldDR,problem, R_num,maxormin,MAXGEN, MAXSIZE, NIND,SUBPOP,GGAP,selectStyle, recombinStyle,recopt,pm, distribute = False)
        elif self.radioButton_q_sorted == True:   #效果较好
            [ObjV, NDSet, NDSetObjV, times] = q_sorted.moea_q_sorted(self,AIM_M, 'aimfuc', PUN_M, 'punishing', FieldDR,problem, R_num,maxormin,MAXGEN, MAXSIZE, NIND,SUBPOP,GGAP,selectStyle, recombinStyle,recopt,pm, distribute = False)
        elif self.radioButton_q_sorted_new == True:   #效果差
            [ObjV, NDSet, NDSetObjV, times] = q_sorted_new.moea_q_sorted_new(self,AIM_M, 'aimfuc', PUN_M, 'punishing', FieldDR,problem, R_num,maxormin,MAXGEN, MAXSIZE, NIND,SUBPOP,GGAP,selectStyle, recombinStyle,recopt,pm, distribute = False)
        else:
            QMessageBox.information(self,'提示信息','请选择算法类型')     
            return
            
        #============================绘制pareto前端===================================
        self.b = NDSet[:, 0:1]
        self.D_z = NDSet[:, 1:2]
        self.d_z = NDSet[:, 2:3]
        self.B = NDSet[:, 3:4]
        self.K_1 = NDSet[:, 4:5]
        self.D_m = NDSet[:, 5:6]
        self.D_r = NDSet[:, 6:7]
        self.L = NDSet[:, 7:8]
        self.z_1 = NDSet[:, 8:9]
        self.z_2 = NDSet[:, 9:10]
        self.m = NDSet[:, 10:11]/4
        self.z_g = NDSet[:, 11:12]*2+1    
        self.Z = NDSet[:, 12:13]
        
        '''
        self.eta = -NDSetObjV[:, 0:1]*100
        self.K_pie = -NDSetObjV[:, 1:2]
        
        plt.rcParams['font.sans-serif'] = ['SimHei']  # 如果要显示中文字体,则在此处设为：SimHei
        
        plt.grid(linestyle="--")  # 设置背景网格线为虚线
        ax = plt.gca()
        ax.spines['top'].set_visible(False)  # 去掉上边框
        ax.spines['right'].set_visible(False)  # 去掉右边框

        plt.scatter(self.eta,self.K_pie)
        
        plt.xlabel('传动效率(%)')
        plt.ylabel('扭转刚度(N·mm/rad)')
        
        plt.show()
        
        
        
        
        
        
        '''
        self.eta = -NDSetObjV[:, 0:1]*100
        self.V = NDSetObjV[:, 1:2]
        self.K_pie = -NDSetObjV[:, 2:]
        
        
        # ===================================绘制三维散点图==============================
         
        # 绘制散点图
        plt.rcParams['font.sans-serif'] = ['SimHei']  # 如果要显示中文字体,则在此处设为：SimHei
        #plt.grid(ls='--')     #设置背景网格线为虚线            
        fig = plt.figure()
        #ax = self.figure.Axes3D(fig)
        ax = Axes3D(fig)
        # 添加坐标轴
        ax.set_xlabel('传动效率(%)')        
        ax.set_ylabel('体积(mm^3)')
        ax.set_zlabel('扭转刚度(N·mm/rad)')

        #ax.grid(False)
        
        # Get rid of colored axes planes
        # First remove fill
        ax.xaxis.pane.fill = False
        ax.yaxis.pane.fill = False
        ax.zaxis.pane.fill = False
        # Now set color to white (or whatever is "invisible")
        ax.xaxis.pane.set_edgecolor('w')
        ax.yaxis.pane.set_edgecolor('w')
        ax.zaxis.pane.set_edgecolor('w')
        #ax.w_xaxis.set_pane_color((1.0, 1.0, 1.0, 1.0))
        ax.scatter(self.eta, self.V, self.K_pie)
        #self.canvas.draw()
        plt.show()

        #ga.frontplot(np.hstack([self.eta, self.V, self.K_pie]), True)
        
        self.tableWidget.setRowCount(NDSet.shape[0])
        for i in range(NDSet.shape[0]):

            item = QTableWidgetItem('%.2f'% (self.b[i:i+1]))
            self.tableWidget.setItem(i, 0, item)            
            item = QTableWidgetItem('%.2f'%(self.D_z[i:i+1, :]))
            self.tableWidget.setItem(i, 1, item)                        
            item = QTableWidgetItem('%.2f'%(self.d_z[i:i+1, :]))
            self.tableWidget.setItem(i, 2, item) 
            item = QTableWidgetItem('%.2f'%(self.B[i:i+1, :]))
            self.tableWidget.setItem(i, 3, item)       
            item = QTableWidgetItem('%.4f'%(self.K_1[i:i+1, :]))
            self.tableWidget.setItem(i, 4, item)        
            item = QTableWidgetItem('%.2f'%(self.D_m[i:i+1, :]))
            self.tableWidget.setItem(i, 5, item)   
            item = QTableWidgetItem('%.2f'%(self.D_r[i:i+1, :]))
            self.tableWidget.setItem(i, 6, item)       
            item = QTableWidgetItem('%.2f'% (self.L[i:i+1, :]))
            self.tableWidget.setItem(i, 7, item)            
            item = QTableWidgetItem('%d'%(self.z_1[i:i+1, :]))
            self.tableWidget.setItem(i, 8, item)                        
            item = QTableWidgetItem('%d'%(self.z_2[i:i+1, :]))
            self.tableWidget.setItem(i, 9, item) 
            item = QTableWidgetItem('%.2f'%(self.m[i:i+1, :]))
            self.tableWidget.setItem(i, 10, item)       
            item = QTableWidgetItem('%d'%(self.z_g[i:i+1, :]))
            self.tableWidget.setItem(i, 11, item)        
            item = QTableWidgetItem('%d'%(self.Z[i:i+1, :]))
            self.tableWidget.setItem(i, 12, item)   
            item = QTableWidgetItem('%.4f'%(self.eta[i:i+1, :]))
            self.tableWidget.setItem(i, 13, item)                 
            item = QTableWidgetItem('%d'%(self.V[i:i+1, :]))
            self.tableWidget.setItem(i, 14, item)       
            item = QTableWidgetItem('%d'%(self.K_pie[i:i+1, :]))
            self.tableWidget.setItem(i, 15, item)     
            
            i_16_5 =  (self.z_g[i:i+1, :]+1)*self.z_2[i:i+1, :]/self.z_1[i:i+1, :]+1
            item = QTableWidgetItem('%.2f'%(i_16_5))
            self.tableWidget.setItem(i, 16, item)        


    @Slot()
    def on_pushButton_select_clicked(self):        
        eta_min = self.doubleSpinBox_eta_min.value()
        V_max = self.doubleSpinBox_V_max.value()
        K_pie_min = self.doubleSpinBox_K_pie_min.value()
        
        self.tableWidget.setRowCount(0)        
        self.tableWidget.setRowCount(len(self.b))
        
        eta = []
        V = []
        K_pie = []      
        j = 0
        for i in range(len(self.b)):
            
            if self.eta[i:i+1, :]>eta_min and self.V[i:i+1, :]<V_max and self.K_pie[i:i+1, :]>K_pie_min:

                item = QTableWidgetItem('%.2f'% (self.b[i:i+1]))
                self.tableWidget.setItem(j, 0, item)            
                item = QTableWidgetItem('%.2f'%(self.D_z[i:i+1, :]))
                self.tableWidget.setItem(j, 1, item)                        
                item = QTableWidgetItem('%.2f'%(self.d_z[i:i+1, :]))
                self.tableWidget.setItem(j, 2, item) 
                item = QTableWidgetItem('%.2f'%(self.B[i:i+1, :]))
                self.tableWidget.setItem(j, 3, item)       
                item = QTableWidgetItem('%.4f'%(self.K_1[i:i+1, :]))
                self.tableWidget.setItem(j, 4, item)        
                item = QTableWidgetItem('%.2f'%(self.D_m[i:i+1, :]))
                self.tableWidget.setItem(j, 5, item)   
                item = QTableWidgetItem('%.2f'%(self.D_r[i:i+1, :]))
                self.tableWidget.setItem(j, 6, item)       
                item = QTableWidgetItem('%.2f'% (self.L[i:i+1, :]))
                self.tableWidget.setItem(j, 7, item)            
                item = QTableWidgetItem('%d'%(self.z_1[i:i+1, :]))
                self.tableWidget.setItem(j, 8, item)                        
                item = QTableWidgetItem('%d'%(self.z_2[i:i+1, :]))
                self.tableWidget.setItem(j, 9, item) 
                item = QTableWidgetItem('%.2f'%(self.m[i:i+1, :]))
                self.tableWidget.setItem(j, 10, item)       
                item = QTableWidgetItem('%d'%(self.z_g[i:i+1, :]))
                self.tableWidget.setItem(j, 11, item)        
                item = QTableWidgetItem('%d'%(self.Z[i:i+1, :]))
                self.tableWidget.setItem(j, 12, item)   
                item = QTableWidgetItem('%.4f'%(self.eta[i:i+1, :]))
                self.tableWidget.setItem(j, 13, item)                 
                item = QTableWidgetItem('%d'%(self.V[i:i+1, :]))
                self.tableWidget.setItem(j, 14, item)       
                item = QTableWidgetItem('%d'%(self.K_pie[i:i+1, :]))
                self.tableWidget.setItem(j, 15, item)        
                
                i_16_5 =  (self.z_g[i:i+1, :]+1)*self.z_2[i:i+1, :]/self.z_1[i:i+1, :]+1
                item = QTableWidgetItem('%.2f'%(i_16_5))
                self.tableWidget.setItem(j, 16, item)    
                
                j = j +1                
                
                eta.append(self.eta[i:i+1, :])
                V.append(self.V[i:i+1, :])
                K_pie.append(self.K_pie[i:i+1, :])
                
        # ===================================绘制三维散点图==============================
         
        # 绘制散点图
        plt.rcParams['font.sans-serif'] = ['SimHei']  # 如果要显示中文字体,则在此处设为：SimHei

        fig = plt.figure()
        #ax = self.figure.Axes3D(fig)
        ax = Axes3D(fig)
        # 添加坐标轴(顺序是Z, Y, X)
        ax.set_zlabel('扭转刚度(N·mm/rad)')
        ax.set_ylabel('体积(mm^3)')
        ax.set_xlabel('传动效率(%)')
        #ax.grid(linestyle = "--")      #设置背景网格线为虚线        
        #.grid(False)
        
        # Get rid of colored axes planes
        # First remove fill
        ax.xaxis.pane.fill = False
        ax.yaxis.pane.fill = False
        ax.zaxis.pane.fill = False
        # Now set color to white (or whatever is "invisible")
        ax.xaxis.pane.set_edgecolor('w')
        ax.yaxis.pane.set_edgecolor('w')
        ax.zaxis.pane.set_edgecolor('w')
        #ax.w_xaxis.set_pane_color((1.0, 1.0, 1.0, 1.0))
        ax.scatter(eta, V, K_pie)
        #self.canvas.draw()
        plt.show()
