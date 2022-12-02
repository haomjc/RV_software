# -*- coding: utf-8 -*-

from PySide2 import QtGui, QtWidgets
from PySide2.QtCore import Slot
from PySide2.QtWidgets import QMainWindow,  QApplication
from Ui_main_window import Ui_MainWindow
from teethmatch_tool_window import Dialog_teethmatch_tool
from GA.GA import Dialog_GA
import sys, math
import numpy as np
import matplotlib.pyplot as plt

class MainWindow(QMainWindow, Ui_MainWindow):

    def __init__(self, parent=None):

        super(MainWindow, self).__init__(parent)
        self.setupUi(self)
        self.Dialog_teethmatch_tool = Dialog_teethmatch_tool()
        self.Dialog_GA = Dialog_GA()
        self.Fw_B = 0.083; self.Fw_C = 1.32; self.Fbeta_A = 0.315; self.Fbeta_C = 1.6; self.fpt_A = 0.063; self.fpt_C = 0.8; self.Fp_A = 0.25; self.Fp_C = 0.63; self.ffalpha_A = 0.063; self.ffalpha_C = 2; self.Fr_A1 = 0.244; self.Fr_C1 = 2.8; self.Fr_A2 = 0.1; self.Fr_C2 = 1.2; self.Fi_A = np.nan; self.Fi_C = np.nan    
        self.radioButton_actual_center_distance = False
        self.radioButton_normal_module = False
        self.speed_ratio_deviation = 0.05    
        
    @Slot()
    def on_action_teeth_match_tool_triggered(self):
        
        self.Dialog_teethmatch_tool.show()
    
    @Slot()
    def on_action_structure_optimal_triggered(self):
        
        self.Dialog_GA.show()
        
    #===============================一级渐开线齿轮 高级===========================================
  
    @Slot()
    def on_pushButton_teethmatch_tool_clicked(self):
                
        self.Dialog_teethmatch_tool.show()
                
    @Slot()
    def on_pushButton_first_calculate_clicked(self):                
                
        first_z1 = self.spinBox_first_z1.value()
        first_z2 = self.spinBox_first_z2.value()
        first_mn = self.doubleSpinBox_first_mn.value()    
        first_ha = self.doubleSpinBox_first_ha.value()   
        first_hf = self.doubleSpinBox_first_hf.value()   
        first_alpha_n = self.doubleSpinBox_first_alpha_n.value()/180*math.pi   
        first_b1 = self.doubleSpinBox_first_b1.value()           
        first_b2 = self.doubleSpinBox_first_b2.value()           
        first_b = self.doubleSpinBox_first_b.value()    
        first_x1 = self.doubleSpinBox_first_x1.value()   
        first_x2 = self.doubleSpinBox_first_x2.value() 
    
        first_d1 = first_mn*first_z1
        first_d2 = first_mn*first_z2
        
        first_c_star = first_hf - first_ha
        
        first_d_b1 = first_d1*math.cos(first_alpha_n)
        first_d_b2 = first_d2*math.cos(first_alpha_n)
        
        h_f1 = (first_ha+first_c_star-first_x1)*first_mn
        h_f2 = (first_ha+first_c_star-first_x2)*first_mn
        
        first_d_f1 = first_d1 - 2*h_f1
        first_d_f2 = first_d2 - 2*h_f2
        
        first_inv_alpha = MainWindow.calc_inv_alpha(self,2*(first_x1+first_x2)*math.tan(first_alpha_n)/(first_z1+first_z2)+math.tan(first_alpha_n)-first_alpha_n)*180/np.pi        
        self.lineEdit_first_inv_alpha.setText('%.2f'%(first_inv_alpha))
        
        first_y = (first_z1+first_z2)/2*(math.cos(first_alpha_n)/math.cos(first_inv_alpha)-1)
        first_a = first_mn*((first_z1+first_z2)/2+first_y)
        
        first_d_a1 = 2*first_a - first_d_f2 - 2*first_c_star*first_mn
        first_d_a2 = 2*first_a - first_d_f1 - 2*first_c_star*first_mn
        
        if first_d_b1/first_d_a1 >=1:
            first_alpha_a_1 = np.arccos(1)
        elif first_d_b1/first_d_a1 <=-1:            
            first_alpha_a_1 = np.arccos(-1)
        else:
            first_alpha_a_1 = np.arccos(first_d_b1/first_d_a1)
            
        if first_d_b2/first_d_a2 >=1:
            first_alpha_a_2 = np.arccos(1)
        elif first_d_b2/first_d_a2 <=-1:            
            first_alpha_a_2 = np.arccos(-1)
        else:
            first_alpha_a_2 = np.arccos(first_d_b2/first_d_a2)
        
        first_epsilon_alpha = 1/(2*math.pi)*(first_z1*(math.tan(first_alpha_a_1)-math.tan(first_inv_alpha))+first_z2*(math.tan(first_alpha_a_2)-math.tan(first_inv_alpha)))
        first_epsilon_beta = 0
        
        first_eta_1 = (math.tan(first_alpha_a_2) - math.tan(first_inv_alpha))/((1+first_z1/first_z2)*math.tan(first_inv_alpha)-math.tan(first_alpha_a_2))*(first_z1+first_z2)/first_z2
        first_eta_2 = (math.tan(first_alpha_a_1) - math.tan(first_inv_alpha))/((1+first_z2/first_z1)*math.tan(first_inv_alpha)-math.tan(first_alpha_a_1))*(first_z1+first_z2)/first_z1        
        
        self.lineEdit_x_sigma.setText(str(first_x1+first_x2))  
        
        self.lineEdit_alpha_a_1.setText('%.2f'%(first_alpha_a_1))        
        self.lineEdit_alpha_a_2.setText('%.2f'%(first_alpha_a_2))
        self.lineEdit_eta_1.setText('%.2f'%(first_eta_1))
        self.lineEdit_eta_2.setText('%.2f'%(first_eta_2))
        self.lineEdit_d_a1.setText('%.2f'%(first_d_a1))
        self.lineEdit_d_a2.setText('%.2f'%(first_d_a2))
        self.lineEdit_d_f1.setText('%.2f'%(first_d_f1))
        self.lineEdit_d_f2.setText('%.2f'%(first_d_f1))
        self.lineEdit_epsilon_alpha.setText('%.2f'%(first_epsilon_alpha))
        self.lineEdit_epsilon_beta.setText('%.2f'%(first_epsilon_beta))    
  
        self.doubleSpinBox_F_w1.setValue(self.Fw_B*math.sqrt(first_d1)+self.Fw_C)
        self.doubleSpinBox_F_w2.setValue(self.Fw_B*math.sqrt(first_d2)+self.Fw_C)  
        
        self.F_beta1 = self.Fbeta_A*math.sqrt(first_b1)+self.Fbeta_C
        self.doubleSpinBox_F_beta1.setValue(self.F_beta1)  
        self.F_beta2 = self.Fbeta_A*math.sqrt(first_b2)+self.Fbeta_C
        self.doubleSpinBox_F_beta2.setValue(self.F_beta2)  
        
        self.fpt_1 = self.fpt_A*first_mn+0.25*self.fpt_A*math.sqrt(first_d1)+self.fpt_C
        self.doubleSpinBox_fpt_1.setValue(self.fpt_1)
        self.fpt_2 = self.fpt_A*first_mn+0.25*self.fpt_A*math.sqrt(first_d2)+self.fpt_C
        self.doubleSpinBox_fpt_2.setValue(self.fpt_2)

        self.doubleSpinBox_Fp_1.setValue(self.Fp_A*math.sqrt(math.pi*first_d1)+self.Fp_C)
        self.doubleSpinBox_Fp_2.setValue(self.Fp_A*math.sqrt(math.pi*first_d2)+self.Fp_C)
        self.doubleSpinBox_ff_alpha_1.setValue(self.ffalpha_A*first_mn+0.0125*self.ffalpha_A*first_d1+self.ffalpha_C)
        self.doubleSpinBox_ff_alpha_2.setValue(self.ffalpha_A*first_mn+0.0125*self.ffalpha_A*first_d2+self.ffalpha_C)
       
        self.doubleSpinBox_first_a.setValue(first_a)   
        
        if self.Fr_A2*first_mn+1.4*self.Fr_A2*math.sqrt(first_d1)+self.Fr_C2 < self.Fr_A1*first_mn+0.25*self.Fr_A1*math.sqrt(first_d1)+self.Fr_C1:
            self.doubleSpinBox_F_r1.setValue(self.Fr_A2*first_mn+1.4*self.Fr_A2*math.sqrt(first_d1)+self.Fr_C2)
        else :
            self.doubleSpinBox_F_r1.setValue(self.Fr_A1*first_mn+0.25*self.Fr_A1*math.sqrt(first_d1)+self.Fr_C1)
        if self.Fr_A2*first_mn+1.4*self.Fr_A2*math.sqrt(first_d2)+self.Fr_C2 < self.Fr_A1*first_mn+0.25*self.Fr_A1*math.sqrt(first_d2)+self.Fr_C1:        
            self.doubleSpinBox_F_r2.setValue(self.Fr_A2*first_mn+1.4*self.Fr_A2*math.sqrt(first_d2)+self.Fr_C2)        
        else:
            self.doubleSpinBox_F_r2.setValue(self.Fr_A1*first_mn+0.25*self.Fr_A1*math.sqrt(first_d2)+self.Fr_C1)
            
        self.lineEdit_F_i1.setText('%.2f'%(self.Fi_A*first_mn+0.25*self.Fi_A*math.sqrt(first_d1)+self.Fi_C))
        self.lineEdit_F_i2.setText('%.2f'%(self.Fi_A*first_mn+0.25*self.Fi_A*math.sqrt(first_d2)+self.Fi_C))
        
    @Slot(str)
    def on_comboBox_accuracy_level_activated(self, p0):

        if p0 == '1级 GB 10095':
            self.Fw_B = 0.083; self.Fw_C = 1.32; self.Fbeta_A = 0.315; self.Fbeta_C = 1.6; self.fpt_A = 0.063; self.fpt_C = 0.8; self.Fp_A = 0.25; self.Fp_C = 0.63; self.ffalpha_A = 0.063; self.ffalpha_C = 2; self.Fr_A1 = 0.244; self.Fr_C1 = 2.8; self.Fr_A2 = 0.1; self.Fr_C2 = 1.2; self.Fi_A = np.nan; self.Fi_C = np.nan
        elif p0 == '2级 GB 10095':
            self.Fw_B = 0.13; self.Fw_C = 2.1; self.Fbeta_A = 0.40; self.Fbeta_C = 2; self.fpt_A = 0.10; self.fpt_C = 1.25; self.Fp_A = 0.4; self.Fp_C = 1; self.ffalpha_A = 0.10; self.ffalpha_C = 2.5; self.Fr_A1 = 0.355; self.Fr_C1 = 4.5; self.Fr_A2 = 0.16; self.Fr_C2 = 1.2; self.Fi_A = np.nan; self.Fi_C = np.nan
        elif p0 == '3级 GB 10095':
            self.Fw_B = 0.21; self.Fw_C = 3.4; self.Fbeta_A = 0.50; self.Fbeta_C = 2.5; self.fpt_A = 0.16; self.fpt_C = 2; self.Fp_A = 0.63; self.Fp_C = 1.6; self.ffalpha_A = 0.16; self.ffalpha_C = 3.15; self.Fr_A1 = 0.56; self.Fr_C1 = 7.1; self.Fr_A2 = 0.25; self.Fr_C2 = 1.2; self.Fi_A = np.nan; self.Fi_C = np.nan    
        elif p0 == '4级 GB 10095':
            self.Fw_B = 0.34; self.Fw_C = 5.4; self.Fbeta_A = 0.63; self.Fbeta_C = 3.15; self.fpt_A = 0.25; self.fpt_C = 3.15; self.Fp_A = 1; self.Fp_C = 2.5; self.ffalpha_A = 0.25; self.ffalpha_C = 4; self.Fr_A1 = 0.90; self.Fr_C1 = 11.2; self.Fr_A2 = 0.4; self.Fr_C2 = 1.2; self.Fi_A = 0.45; self.Fi_C = 5.6            
        elif p0 == '5级 GB 10095': 
            self.Fw_B = 0.54; self.Fw_C = 8.7; self.Fbeta_A = 0.8; self.Fbeta_C = 4; self.fpt_A = 0.40; self.fpt_C = 5; self.Fp_A = 1.6; self.Fp_C = 4; self.ffalpha_A = 0.4; self.ffalpha_C = 5; self.Fr_A1 = 1.40; self.Fr_C1 = 18; self.Fr_A2 = 0.63; self.Fr_C2 = 1.2; self.Fi_A = 0.63; self.Fi_C = 8          
        elif p0 == '6级 GB 10095':  
            self.Fw_B = 0.87; self.Fw_C = 14; self.Fbeta_A = 1; self.Fbeta_C = 5; self.fpt_A = 0.63; self.fpt_C = 8; self.Fp_A = 2.5; self.Fp_C = 6.3; self.ffalpha_A = 0.63; self.ffalpha_C = 6.3; self.Fr_A1 = 2.24; self.Fr_C1 = 28; self.Fr_A2 = 1; self.Fr_C2 = 1.2; self.Fi_A = 0.9; self.Fi_C = 11.2           
        elif p0 == '7级 GB 10095':    
            self.Fw_B = 1.22; self.Fw_C = 19.4; self.Fbeta_A = 1.25; self.Fbeta_C = 6.3; self.fpt_A = 0.90; self.fpt_C = 11.2; self.Fp_A = 3.55; self.Fp_C = 9; self.ffalpha_A = 1; self.ffalpha_C = 8; self.Fr_A1 = 3.15; self.Fr_C1 = 40; self.Fr_A2 = 1.4; self.Fr_C2 = 1.2; self.Fi_A = 1.25; self.Fi_C = 16           
        elif p0 == '8级 GB 10095':   
            self.Fw_B = 1.7; self.Fw_C = 27; self.Fbeta_A = 2; self.Fbeta_C = 10; self.fpt_A = 1.25; self.fpt_C = 16; self.Fp_A = 5; self.Fp_C = 12.5; self.ffalpha_A = 1.6; self.ffalpha_C = 10; self.Fr_A1 = 4; self.Fr_C1 = 50; self.Fr_A2 = 1.75; self.Fr_C2 = 1.2; self.Fi_A = 1.8; self.Fi_C = 22.4          
        elif p0 == '9级 GB 10095':
            self.Fw_B = 2.4; self.Fw_C = 38; self.Fbeta_A = 3.15; self.Fbeta_C = 16; self.fpt_A = 1.8; self.fpt_C = 22.4; self.Fp_A = 7.1; self.Fp_C = 18; self.ffalpha_A = 2.5; self.ffalpha_C = 16; self.Fr_A1 = 5; self.Fr_C1 = 63; self.Fr_A2 = np.nan; self.Fr_C2 = np.nan; self.Fi_A = 2.24; self.Fi_C = 28         
        elif p0 == '10级 GB 10095':
            self.Fw_B = 3.3; self.Fw_C = 53; self.Fbeta_A = 5; self.Fbeta_C = 25; self.fpt_A = 2.5; self.fpt_C = 31.5; self.Fp_A = 10; self.Fp_C = 25; self.ffalpha_A = 4; self.ffalpha_C = 25; self.Fr_A1 = 6.3; self.Fr_C1 = 80; self.Fr_A2 = np.nan; self.Fr_C2 = np.nan; self.Fi_A = 2.8; self.Fi_C = 35.5         
        elif p0 == '11级 GB 10095':
            self.Fw_B = 4.7; self.Fw_C = 74; self.Fbeta_A = 8; self.Fbeta_C = 40; self.fpt_A = 3.55; self.fpt_C = 45; self.Fp_A = 14; self.Fp_C = 35.5; self.ffalpha_A = 6.3; self.ffalpha_C = 40; self.Fr_A1 = 8; self.Fr_C1 = 100; self.Fr_A2 = np.nan; self.Fr_C2 = np.nan; self.Fi_A = 3.55; self.Fi_C = 45         
        elif p0 == '12级 GB 10095':
            self.Fw_B = 6.5; self.Fw_C = 104; self.Fbeta_A = 12.5; self.Fbeta_C = 63; self.fpt_A = 5; self.fpt_C = 63; self.Fp_A = 20; self.Fp_C = 50; self.ffalpha_A = 10; self.ffalpha_C = 63; self.Fr_A1 = 10; self.Fr_C1 = 125; self.Fr_A2 = np.nan; self.Fr_C2 = np.nan; self.Fi_A = 4.5; self.Fi_C = 56


    def calc_inv_alpha(self, inv):
        alpha_low = 0
        alpha_high = 2
        
        for i in range(30):
            alpha_mid = (alpha_low+alpha_high)/2
            if math.tan(alpha_mid) - alpha_mid > inv :
                alpha_high = alpha_mid
            else :
                alpha_low = alpha_mid
                
        return alpha_mid
        
    #================================================一级渐开线齿轮 校核======================================================      
    @Slot()
    def on_pushButton_first_check_clicked(self): 
        
        first_z1 = self.spinBox_first_z1.value()
        first_z2 = self.spinBox_first_z2.value()
        first_mn = self.doubleSpinBox_first_mn.value()    
        first_ha = self.doubleSpinBox_first_ha.value()   
        first_hf = self.doubleSpinBox_first_hf.value()   
        first_alpha_n = self.doubleSpinBox_first_alpha_n.value()/180*math.pi   
        first_b1 = self.doubleSpinBox_first_b1.value()           
        first_b2 = self.doubleSpinBox_first_b2.value()           
        first_b = self.doubleSpinBox_first_b.value()    
        first_x1 = self.doubleSpinBox_first_x1.value()   
        first_x2 = self.doubleSpinBox_first_x2.value() 
        first_rpm = self.doubleSpinBox_first_rpm.value()
        first_torsion = self.doubleSpinBox_first_torsion.value()
        first_planetary_wheel = self.spinBox_planetary_wheel.value()
        first_ka = self.doubleSpinBox_first_ka.value()
        delta_hlim1 = self.doubleSpinBox_delta_hlim1.value()
        delta_hlim2 = self.doubleSpinBox_delta_hlim2.value()
        delta_flim1 = self.doubleSpinBox_delta_flim1.value()
        delta_flim2 = self.doubleSpinBox_delta_flim2.value()
        first_v40 = self.spinBox_first_v40.value()
        first_k = self.doubleSpinBox_first_k.value()
        first_l = self.doubleSpinBox_first_l.value()
        first_s = self.doubleSpinBox_first_s.value()
        first_dsh = self.doubleSpinBox_first_dsh.value()
        first_e1 = self.doubleSpinBox_e1.value()
        first_e2 = self.doubleSpinBox_e2.value()  
        first_nu1 = self.doubleSpinBox_nu1.value()        
        first_nu2 = self.doubleSpinBox_nu2.value()
        hb1 = self.doubleSpinBox_hb1.value()
        hb2 = self.doubleSpinBox_hb2.value()
        z_nt1 = self.doubleSpinBox_z_nt1.value()
        z_nt2 = self.doubleSpinBox_z_nt2.value()
        y_nt1 = self.doubleSpinBox_y_nt1.value()
        y_nt2 = self.doubleSpinBox_y_nt2.value()        
        s_hmin1 = self.doubleSpinBox_s_hmin1.value()
        s_hmin2 = self.doubleSpinBox_s_hmin2.value()
        s_fmin1 = self.doubleSpinBox_s_fmin1.value()
        s_fmin2 = self.doubleSpinBox_s_fmin2.value()        
        first_rho_fp = self.doubleSpinBox_first_rho_fp.value()
        rho_pie1 = self.doubleSpinBox_rho_pie1.value()
        rho_pie2 = self.doubleSpinBox_rho_pie2.value()
        
        first_d1 = first_mn*first_z1
        first_d2 = first_mn*first_z2
        
        first_c_star = first_hf - first_ha
        
        first_d_b1 = first_d1*math.cos(first_alpha_n)
        first_d_b2 = first_d2*math.cos(first_alpha_n)
        
        h_f1 = (first_ha+first_c_star-first_x1)*first_mn
        h_f2 = (first_ha+first_c_star-first_x2)*first_mn
        
        first_d_f1 = first_d1 - 2*h_f1
        first_d_f2 = first_d2 - 2*h_f2
        
        first_inv_alpha = MainWindow.calc_inv_alpha(self,2*(first_x1+first_x2)*math.tan(first_alpha_n)/(first_z1+first_z2)+math.tan(first_alpha_n)-first_alpha_n)*180/np.pi        
        self.lineEdit_first_inv_alpha.setText('%.2f'%(first_inv_alpha))
        
        first_y = (first_z1+first_z2)/2*(math.cos(first_alpha_n)/math.cos(first_inv_alpha)-1)
        first_a = first_mn*((first_z1+first_z2)/2+first_y)
        
        first_d_a1 = 2*first_a - first_d_f2 - 2*first_c_star*first_mn
        first_d_a2 = 2*first_a - first_d_f1 - 2*first_c_star*first_mn
        
        if first_d_b1/first_d_a1 >=1:
            first_alpha_a_1 = np.arccos(1)
        elif first_d_b1/first_d_a1 <=-1:            
            first_alpha_a_1 = np.arccos(-1)
        else:
            first_alpha_a_1 = np.arccos(first_d_b1/first_d_a1)
            
        if first_d_b2/first_d_a2 >=1:
            first_alpha_a_2 = np.arccos(1)
        elif first_d_b2/first_d_a2 <=-1:            
            first_alpha_a_2 = np.arccos(-1)
        else:
            first_alpha_a_2 = np.arccos(first_d_b2/first_d_a2)
        
        first_epsilon_alpha = 1/(2*math.pi)*(first_z1*(math.tan(first_alpha_a_1)-math.tan(first_inv_alpha))+first_z2*(math.tan(first_alpha_a_2)-math.tan(first_inv_alpha)))
        first_epsilon_beta = 0
        
        first_eta_1 = (math.tan(first_alpha_a_2) - math.tan(first_inv_alpha))/((1+first_z1/first_z2)*math.tan(first_inv_alpha)-math.tan(first_alpha_a_2))*(first_z1+first_z2)/first_z2
        first_eta_2 = (math.tan(first_alpha_a_1) - math.tan(first_inv_alpha))/((1+first_z2/first_z1)*math.tan(first_inv_alpha)-math.tan(first_alpha_a_1))*(first_z1+first_z2)/first_z1        
        
        #==========================================齿面接触强度核算====================================
        
        first_F_t = 2000*first_torsion/(first_planetary_wheel*first_mn*first_z1)
        
        first_c1 = -0.5048*math.log(first_z1) -1.144*math.log(first_mn) + 2.852*math.log(self.fpt_1) + 3.32
        first_c2 = -0.5048*math.log(first_z2) -1.144*math.log(first_mn) + 2.852*math.log(self.fpt_2) + 3.32
        
        first_b1 = 0.25*(first_c1-5.0)**0.667
        first_b2 = 0.25*(first_c2-5.0)**0.667
        
        first_a1 = 50 +56*(1.0-first_b1)    
        first_a2 = 50 +56*(1.0-first_b2) 
        
        first_v1 = first_d1*math.pi*first_rpm/60
        first_kv1 = (first_a1/(first_a1+math.sqrt(200*first_v1)))**(-first_b1)
        first_kv2 = (first_a2/(first_a2+math.sqrt(200*first_v1)))**(-first_b2)    
      
        first_wm1 = first_F_t*first_ka*first_kv1/first_b
        first_wm2 = first_F_t*first_ka*first_kv2/first_b
        
        first_x_beta1 = 1 - 320/delta_hlim1
        first_x_beta2 = 1 - 320/delta_hlim2
        
        first_gamma1 = (abs(1 + first_k*first_l*first_s/(first_d1**2)*(first_d1/first_dsh)**4-0.3)+0.3)*(first_b/first_d1)**2
        first_gamma2 = (abs(1 + first_k*first_l*first_s/(first_d2**2)*(first_d2/first_dsh)**4-0.3)+0.3)*(first_b/first_d2)**2      
        first_f_sh1 = first_wm1*0.023*first_gamma1
        first_f_sh2 = first_wm2*0.023*first_gamma2   
        first_F_beta_x1 = abs(1.33*first_f_sh1 - self.F_beta1)
        first_F_beta_x2 = abs(1.33*first_f_sh2 - self.F_beta2)
        first_F_beta_y1 = first_F_beta_x1*first_x_beta1
        first_F_beta_y2 = first_F_beta_x2*first_x_beta2

        first_c_gamma = 20
        
        if math.sqrt(2*first_wm1/(first_F_beta_y1*first_c_gamma)) <= 1:
            first_k_h_beta1 = math.sqrt(2*first_F_beta_y1*first_c_gamma/first_wm1)
        else:
            first_k_h_beta1 = 1+0.5*first_F_beta_y1*first_c_gamma/first_wm1
        if math.sqrt(2*first_wm2/(first_F_beta_y2*first_c_gamma)) <= 1:
            first_k_h_beta2 = math.sqrt(2*first_F_beta_y2*first_c_gamma/first_wm2)
        else:
            first_k_h_beta2 = 1+0.5*first_F_beta_y2*first_c_gamma/first_wm2
            
        first_f_th1 = first_F_t*first_ka*first_kv1*first_k_h_beta1
        first_f_th2 = first_F_t*first_ka*first_kv2*first_k_h_beta2

        first_epsilon_gamma = first_epsilon_alpha+first_epsilon_beta
        
        first_y_alpha = 160/delta_hlim1*self.fpt

        if first_epsilon_gamma <= 2:
            first_k_h_alpha1 = first_epsilon_gamma/2*(0.9+0.4*first_c_gamma*(self.fpt_1-first_y_alpha)/(first_f_th1/first_b))
            first_k_h_alpha2 = first_epsilon_gamma/2*(0.9+0.4*first_c_gamma*(self.fpt_2-first_y_alpha)/(first_f_th2/first_b))
        else:
            first_k_h_alpha1 = 0.9+0.4*math.sqrt(2*(first_epsilon_gamma-1)/first_epsilon_gamma)*first_c_gamma*(self.fpt_1-first_y_alpha)/(first_f_th1/first_b)
            first_k_h_alpha2 = 0.9+0.4*math.sqrt(2*(first_epsilon_gamma-1)/first_epsilon_gamma)*first_c_gamma*(self.fpt_2-first_y_alpha)/(first_f_th2/first_b)          
          
        if first_k_h_alpha1<1.0:
            first_k_h_alpha1 = 1.0
            
        if first_k_h_alpha2 < 1.0:
            first_k_h_alpha2 = 1.0

        first_m1 = math.tan(first_inv_alpha)/math.sqrt((math.sqrt(first_d_a1**2/(first_d_b1)-1)-2*math.pi*first_z1)*(math.sqrt(first_d_a2**2/(first_d_b2)-1)-(first_epsilon_alpha-1)*2*math.pi/first_z2))
        first_m2 = math.tan(first_inv_alpha)/math.sqrt((math.sqrt(first_d_a2**2/(first_d_b2)-1)-2*math.pi*first_z2)*(math.sqrt(first_d_a1**2/(first_d_b1)-1)-(first_epsilon_alpha-1)*2*math.pi/first_z1))        
        if first_m1 >1:
            first_zb = first_m1
        else:
            first_zb = 1
        if first_m2 >1:
            first_zd = first_m2
        else:
            first_zd = 1
            
        first_zh =  math.sqrt(2*math.cos(first_inv_alpha)/((math.cos(first_alpha_n)**2)*math.sin(first_inv_alpha)))           
        first_ze = math.sqrt(1/(math.pi*((1-first_nu1)/first_e1+(1-first_nu2)/first_e2)))
        first_z_epsilon = math.sqrt((4-first_epsilon_alpha)/3)
        first_z_beta = 1
        first_delta_h0 = first_zh*first_ze*first_z_epsilon*first_z_beta*math.sqrt(first_F_t*(first_z1+first_z2)/(first_d1*first_b*first_z2))
        first_delta_h1 = first_zb*first_delta_h0*math.sqrt(first_ka*first_kv1*first_k_h_beta1*first_k_h_alpha1)
        first_delta_h2 = first_zd*first_delta_h0*math.sqrt(first_ka*first_kv2*first_k_h_beta2*first_k_h_alpha2)     
        
        z_l_z_v_z_r = 1
        z_x = 1
        z_w1= 1.2-(hb1-130)/1700
        z_w2= 1.2-(hb2-130)/1700        
        
        first_delta_hg1 = delta_hlim1*z_nt1*z_l_z_v_z_r*z_w1*z_x
        first_delta_hg2 = delta_hlim2*z_nt2*z_l_z_v_z_r*z_w2*z_x   
   
        first_delta_hp1 = first_delta_hg1/s_hmin1
        first_delta_hp2 = first_delta_hg2/s_hmin2
        first_sh_1 = first_delta_hg1/first_delta_h1
        first_sh_2 = first_delta_hg2/first_delta_h2       
        
        self.lineEdit_delta_h1.setText('%.2f'%(first_delta_h1))
        self.lineEdit_delta_h2.setText('%.2f'%(first_delta_h2))
        self.lineEdit_delta_hp1.setText('%.2f'%(first_delta_hp1))        
        self.lineEdit_delta_hp2.setText('%.2f'%(first_delta_hp2))        
        self.lineEdit_sh_1.setText('%.2f'%(first_sh_1))        
        self.lineEdit_sh_2.setText('%.2f'%(first_sh_2))        
        
        #============================================轮齿弯曲强度核算===========================================
        
        first_hfp = first_mn*first_hf
        first_s_pr = 0
        first_E_1 = math.pi*first_mn/4-first_hfp*math.tan(first_alpha_n)+first_s_pr/math.cos(first_alpha_n)-(1-math.sin(first_alpha_n))*first_rho_fp/math.cos(first_alpha_n)
        first_G_1 = first_rho_fp/first_mn-first_hfp/first_mn+first_x1
        first_H_1 = 2/first_z1*(math.pi/2-first_E_1/first_mn)-math.pi/3
        first_theta_1 = MainWindow.calc_first_theta(self, first_G_1, first_z1, first_H_1)
        first_s_fnmn1 = first_z1*math.sin(math.pi/3 - first_theta_1)+math.sqrt(3)*(first_G_1/math.cos(first_theta_1)-first_rho_fp/first_mn)
        first_rho_fmn1 = first_rho_fp/first_mn+2*first_G_1**2/(math.cos(first_theta_1)*(first_z1*math.cos(first_theta_1)**2-2*first_G_1))
        first_epsilon_alpha_n = first_epsilon_alpha
        first_d_n1 = first_mn*first_z1
        first_d_bn1 = first_d_n1*math.cos(first_alpha_n)
        first_d_an1 = first_d_n1+first_d_a1-first_d1
        first_d_en1 = 2*math.sqrt((math.sqrt((first_d_an1/2)**2+(first_d_bn1/2)**2)-math.pi*first_mn*math.cos(first_alpha_n)*(first_epsilon_alpha_n-1))**2+(first_d_bn1/2)**2)
        first_alpha_en1 = math.acos(first_d_bn1/first_d_en1)
        first_gamma_e1 = 1/first_z1*(math.pi/2+2*first_x1*math.tan(first_alpha_n))+first_alpha_n-math.tan(first_alpha_n)-(first_alpha_en1-math.tan(first_alpha_en1))
        first_alpha_fen1 = first_alpha_en1-first_gamma_e1
        first_h_femn1 = 1/2*(math.cos(first_gamma_e1)-math.sin(first_gamma_e1)*math.tan(first_alpha_fen1)*first_d_en1/first_mn-first_z1*math.cos(math.pi/3-first_theta_1)-first_G_1/math.cos(first_theta_1)+first_rho_fp/first_mn)
        first_Y_F1 = 6*first_h_femn1*math.cos(first_alpha_fen1)/(first_s_fnmn1**2*math.cos(first_alpha_n))
        first_L1 = first_s_fnmn1/first_h_femn1
        first_q_s1 = first_s_fnmn1/first_rho_fmn1        
        first_Y_S1 = (1.2+0.13*first_L1)*(1/2*first_q_s1)**(1/(1.121+2.3/first_L1))        
        
        first_hfp = first_mn*first_hf
        first_s_pr = 0
        first_E_2 = math.pi*first_mn/4-first_hfp*math.tan(first_alpha_n)+first_s_pr/math.cos(first_alpha_n)-(1-math.sin(first_alpha_n))*first_rho_fp/math.cos(first_alpha_n)
        first_G_2 = first_rho_fp/first_mn-first_hfp/first_mn+first_x2
        first_H_2 = 2/first_z1*(math.pi/2-first_E_2/first_mn)-math.pi/3
        first_theta_2 = MainWindow.calc_first_theta(self, first_G_2, first_z2, first_H_2)
        first_s_fnmn2 = first_z2*math.sin(math.pi/3 - first_theta_2)+math.sqrt(3)*(first_G_2/math.cos(first_theta_2)-first_rho_fp/first_mn)
        first_rho_fmn2 = first_rho_fp/first_mn+2*first_G_2**2/(math.cos(first_theta_2)*(first_z2*math.cos(first_theta_2)**2-2*first_G_2))
        first_epsilon_alpha_n = first_epsilon_alpha
        first_d_n2 = first_mn*first_z2
        first_d_bn2 = first_d_n2*math.cos(first_alpha_n)
        first_d_an2 = first_d_n2+first_d_a2-first_d2
        first_d_en2 = 2*math.sqrt((math.sqrt((first_d_an2/2)**2+(first_d_bn2/2)**2)-math.pi*first_mn*math.cos(first_alpha_n)*(first_epsilon_alpha_n-1))**2+(first_d_bn2/2)**2)
        first_alpha_en2 = math.acos(first_d_bn2/first_d_en2)
        first_gamma_e2 = 1/first_z2*(math.pi/2+2*first_x2*math.tan(first_alpha_n))+first_alpha_n-math.tan(first_alpha_n)-(first_alpha_en2-math.tan(first_alpha_en2))
        first_alpha_fen2 = first_alpha_en2-first_gamma_e2
        first_h_femn2 = 1/2*(math.cos(first_gamma_e2)-math.sin(first_gamma_e2)*math.tan(first_alpha_fen2)*first_d_en2/first_mn-first_z2*math.cos(math.pi/3-first_theta_2)-first_G_2/math.cos(first_theta_2)+first_rho_fp/first_mn)
        first_Y_F2 = 6*first_h_femn2*math.cos(first_alpha_fen2)/(first_s_fnmn2**2*math.cos(first_alpha_n))  
        first_L2 = first_s_fnmn2/first_h_femn2
        first_q_s2 = first_s_fnmn2/first_rho_fmn2
        first_Y_S2 = (1.2+0.13*first_L2)*(1/2*first_q_s2)**(1/(1.121+2.3/first_L2))     
        
        first_Y_beta = 1
        
        first_delta_F0_1 = first_F_t/(first_b*first_mn)*first_Y_F1*first_Y_S1*first_Y_beta
        first_delta_F0_2 = first_F_t/(first_b*first_mn)*first_Y_F2*first_Y_S2*first_Y_beta     
     
        first_h = first_ha+first_hf
        
        if first_b1> first_b2:
            first_N = (first_b2/first_h)**2/(1 +(first_b2/first_h)+(first_b2/first_h)**2)
        else:
            first_N = (first_b1/first_h)**2/(1 +(first_b1/first_h)+(first_b1/first_h)**2)
        first_k_f_beta1 = first_k_h_beta1**first_N
        first_k_f_beta2 = first_k_h_beta2**first_N 
        first_k_f_alpha1 = first_k_h_alpha1
        first_k_f_alpha2 = first_k_h_alpha2
        
        first_delta_f1 = first_delta_F0_1*first_ka*first_kv1*first_k_f_beta1*first_k_f_alpha1
        first_delta_f2 = first_delta_F0_2*first_ka*first_kv2*first_k_f_beta2*first_k_f_alpha2
        
        first_Y_ST = 2.0
        first_q_sT = 2.5
        first_X_star_1 = 1/5*(1+2*first_q_s1)
        first_X_star_T = 1/5*(1+2*first_q_sT)
        
        first_Y_delta_rel_t1 = (1+math.sqrt(rho_pie1*first_X_star_1))/(1+math.sqrt(rho_pie1*first_X_star_T))
        first_Rz = 3.2
        first_Y_R_rel_t = 1.674-0.529*(first_Rz+1)**0.1
        first_Y_x = 1.0
        first_delta_fg1 = delta_flim1*first_Y_ST*y_nt1*first_Y_delta_rel_t1*first_Y_R_rel_t*first_Y_x
        first_delta_fp1 = first_delta_fg1/s_fmin1
        first_sf_1 = first_delta_fg1/first_delta_f1
        
        first_X_star_2 = 1/5*(1+2*first_q_s2)
        first_Y_delta_rel_t2 = (1+math.sqrt(rho_pie2*first_X_star_2))/(1+math.sqrt(rho_pie2*first_X_star_T))
        first_delta_fg2 = delta_flim2*first_Y_ST*y_nt2*first_Y_delta_rel_t2*first_Y_R_rel_t*first_Y_x
        first_delta_fp2 = first_delta_fg2/s_fmin2
        first_sf_2 = first_delta_fg2/first_delta_f2    
       
        self.lineEdit_delta_f1.setText('%.2f'%(first_delta_f1))
        self.lineEdit_delta_f2.setText('%.2f'%(first_delta_f2))
        self.lineEdit_delta_fp1.setText('%.2f'%(first_delta_fp1))        
        self.lineEdit_delta_fp2.setText('%.2f'%(first_delta_fp2))        
        self.lineEdit_sf_1.setText('%.2f'%(first_sf_1))        
        self.lineEdit_sf_2.setText('%.2f'%(first_sf_2))            
        
    def calc_first_theta(self, G, Zn, H):
        theta = -H/(1-2*G/Zn)
        for i in range(20):
            theta = 2*G/Zn*math.tan(theta)-H
        return theta
    
    #=====================================二级摆线针轮 高级================================================    
    @Slot()
    def on_pushButton_second_calculate_clicked(self):
        
        second_zp = self.spinBox_second_zp.value()
        second_rp = self.doubleSpinBox_second_rp.value()
        second_k1_origin = self.doubleSpinBox_second_k1_origin.value()        
        second_k2_origin = self.doubleSpinBox_second_k2_origin.value()     
        second_b = self.doubleSpinBox_second_b.value()            

        second_zc = second_zp - 1
        second_a = second_k1_origin*second_rp/second_zp
        second_r_rp = second_rp/second_k2_origin*math.sin(math.pi/second_zp)
        
        self.doubleSpinBox_second_a.setValue(second_a)
        self.doubleSpinBox_second_r_rp.setValue(second_r_rp)
        self.spinBox_second_zc.setValue(second_zc)    
        
    @Slot()
    def on_pushButton_second_set_clicked(self):
        second_zp = self.spinBox_second_zp.value()
        second_rp = self.doubleSpinBox_second_rp.value()
        second_a = self.doubleSpinBox_second_a.value()     
        second_r_rp = self.doubleSpinBox_second_r_rp.value()

        second_k1 = second_a*second_zp/second_rp
        second_k2 = second_rp/second_r_rp*math.sin(math.pi/second_zp)
        second_zc = second_zp - 1
        second_a = second_k1*second_rp/second_zp
        second_rp_pie = second_a*second_zp
        second_rc_pie = second_a*second_zc
        second_rp = second_a*second_zp/second_k1
        second_rc = second_a*second_zc/second_k1
        second_r_ap = second_rp - second_r_rp
        second_r_fp = second_rp
        second_r_ac = second_rp + second_a - second_r_rp
        second_r_fc = second_rp - second_a - second_r_rp        
        second_p_pie = 2*math.pi*second_a
        
        self.lineEdit_second_k1.setText('%.2f'%(second_k1))
        self.lineEdit_second_k2.setText('%.2f'%(second_k2))
        self.lineEdit_second_rp_pie.setText('%.2f'%(second_rp_pie))
        self.lineEdit_second_rc_pie.setText('%.2f'%(second_rc_pie))
        self.lineEdit_second_rc.setText('%.2f'%(second_rc))
        self.lineEdit_second_r_ap.setText('%.2f'%(second_r_ap))
        self.lineEdit_second_r_fp.setText('%.2f'%(second_r_fp))
        self.lineEdit_second_r_ac.setText('%.2f'%(second_r_ac))
        self.lineEdit_second_r_fc.setText('%.2f'%(second_r_fc))
        self.lineEdit_second_p_pie.setText('%.2f'%(second_p_pie))
        self.spinBox_second_zc.setValue(second_zc)    
        
    #===========================================绘制摆线轮曲线================================================
    @Slot()
    def on_pushButton_second_draw_cycloid_gear_clicked(self):
        second_zp = self.spinBox_second_zp.value()
        second_rp = self.doubleSpinBox_second_rp.value()
        second_a = self.doubleSpinBox_second_a.value()  
        second_r_rp = self.doubleSpinBox_second_r_rp.value()
        second_delta_rp = self.doubleSpinBox_second_delta_rp.value()
        second_delta_r_rp = self.doubleSpinBox_second_delta_r_rp.value()
        second_delta = self.doubleSpinBox_second_delta.value()/180*math.pi

        second_phi_b = np.arange(0, 2*(second_zp-1)*math.pi, 0.0001)        
        second_k1_pie = second_a*second_zp/(second_rp - second_delta_rp)
        second_s_pie = 1+second_k1_pie**2-2*second_k1_pie*np.cos(second_phi_b)
        second_i_H = second_zp/(second_zp-1)
        second_x = -((second_rp - second_delta_rp)-(second_r_rp+second_delta_r_rp)*second_s_pie**(-1/2))*np.sin((1-second_i_H)*second_phi_b-second_delta) - second_a/(second_rp - second_delta_rp)*((second_rp - second_delta_rp)-second_zp*(second_r_rp+second_delta_r_rp)*second_s_pie**(-1/2))*np.sin(second_i_H*second_phi_b+second_delta)
        second_y = ((second_rp - second_delta_rp)-(second_r_rp+second_delta_r_rp)*second_s_pie**(-1/2))*np.cos((1-second_i_H)*second_phi_b-second_delta) - second_a/(second_rp - second_delta_rp)*((second_rp - second_delta_rp)-second_zp*(second_r_rp+second_delta_r_rp)*second_s_pie**(-1/2))*np.cos(second_i_H*second_phi_b+second_delta)   
        '''
        second_cos_beta = (second_k1_pie*np.sin(second_zp*second_phi_b) - np.sin(second_phi_b))/np.sqrt(1+second_k1_pie**2-2*second_k1_pie*np.cos(second_zp*second_phi_b))
        second_sin_beta = (-second_k1_pie*np.cos(second_zp*second_phi_b) + np.cos(second_phi_b))/np.sqrt(1+second_k1_pie**2-2*second_k1_pie*np.cos(second_zp*second_phi_b))
        second_x = (second_rp - second_delta_rp)*(np.sin(second_phi_b+second_delta) - (second_k1_pie/second_zp)*np.sin(second_zp*second_phi_b+second_delta)) + (second_r_rp+second_delta_r_rp)*(second_cos_beta*np.cos(second_delta)-second_sin_beta*np.sin(second_delta))
        second_y = (second_rp - second_delta_rp)*(np.cos(second_phi_b+second_delta) - (second_k1_pie/second_zp)*np.cos(second_zp*second_phi_b+second_delta)) - (second_r_rp+second_delta_r_rp)*(second_sin_beta*np.cos(second_delta)+second_cos_beta*np.sin(second_delta)) 
        '''
        ax = plt.subplot(111)
        ax.set_aspect(1)
        ax.plot(second_x, second_y)
        plt.show()
        
    #总体设计    
    @Slot()
    def on_pushButton_main_zp_calc_clicked(self): 
        gearteeth_min = self.spinBox_gearteeth_min.value()
        first_speed_ratio = self.doubleSpinBox_first_speed_ratio.value()
        main_speed_ratio = self.doubleDoubleSpinBox_main_speed_ratio.value()
        
        main_zp = round((main_speed_ratio-1)/first_speed_ratio/2, 0)*2
        self.spinBox_main_zp.setValue(main_zp)
        
    @Slot()
    def on_pushButton_main_a0_calc_clicked(self):    
        main_rp = self.doubleSpinBox_main_rp.value()
        
        main_a0_min = str(0.5*main_rp)
        main_a0_max = str(0.6*main_rp)
        self.label_main_a0.setText(QtWidgets.QApplication.translate("MainWindow", main_a0_min+"~"+main_a0_max, None, -1))
        
    @Slot()
    def on_pushButton_main_speed_ratio_calc_clicked(self): 
        main_zp =  self.spinBox_main_zp.value()
        main_z1 = self.spinBox_main_z1.value()
        main_z2 = self.spinBox_main_z2.value()

        main_real_speed_ratio = 1+main_z2/main_z1*main_zp
        self.doubleSpinBox_main_real_speed_ratio.setValue(main_real_speed_ratio)


    @Slot(str)
    def on_comboBox_speed_ratio_deviation_activated(self, p0):

        if p0 == '±5%':
            self.speed_ratio_deviation = 0.05
        elif p0 == '±10%':
            self.speed_ratio_deviation = 0.10    
        elif p0 == '±15%':
            self.speed_ratio_deviation = 0.15        
        elif p0 == '±20%':
            self.speed_ratio_deviation = 0.20   
            
    @Slot()
    def on_radioButton_actual_center_distance_clicked(self):
        self.radioButton_actual_center_distance = True
        self.radioButton_normal_module = False     
        
    @Slot()
    def on_radioButton_normal_module_clicked(self):
        self.radioButton_normal_module = True    
        self.radioButton_actual_center_distance = False
        
    @Slot()
    def on_pushButton_main_calculate_clicked(self):    
        self.speed_ratio = self.doubleSpinBox_first_speed_ratio.value()
        self.gearteeth_min = self.spinBox_gearteeth_min.value()
        self.d_or_mn = self.doubleSpinBox_d_or_mn.value()
        self.main_planetary_wheel = self.spinBox_main_planetary_wheel.value()
        
        if self.radioButton_actual_center_distance == True:
            self.actual_center_distance = self.d_or_mn
        elif self.radioButton_normal_module == True:
            self.normal_module = self.d_or_mn
            
        i = 0
        self.tableWidget.setRowCount(int((1+self.speed_ratio_deviation)*self.gearteeth_min*self.speed_ratio)-math.ceil((1-self.speed_ratio_deviation)*self.gearteeth_min*self.speed_ratio))

        for gearteeth_max in range(math.ceil((1-self.speed_ratio_deviation)*self.gearteeth_min*self.speed_ratio),int((1+self.speed_ratio_deviation)*self.gearteeth_min*self.speed_ratio) ):

            if (self.main_planetary_wheel  == 2 and (gearteeth_max - self.gearteeth_min)%2 ==0) or (self.main_planetary_wheel  == 3 and gearteeth_max%2 ==0 and (2*self.gearteeth_min - gearteeth_max)%6 == 0) or (self.main_planetary_wheel  == 3 and gearteeth_max%2 ==1 and (2*self.gearteeth_min - gearteeth_max )%6 == 3):

                speed_ratio_deviation = gearteeth_max/self.gearteeth_min
                if self.radioButton_actual_center_distance == True:
                    self.normal_module = 2*self.actual_center_distance/(self.gearteeth_min*(speed_ratio_deviation+1))
                elif self.radioButton_normal_module == True:
                    self.actual_center_distance = self.gearteeth_min*self.normal_module*(speed_ratio_deviation+1)/2
    
                item = QtWidgets.QTableWidgetItem('%d'% (self.gearteeth_min))
                self.tableWidget.setItem(i, 0, item)            
                item = QtWidgets.QTableWidgetItem('%d'%(gearteeth_max))
                self.tableWidget.setItem(i, 1, item)                        
                item = QtWidgets.QTableWidgetItem('%.2f'%(self.normal_module))
                self.tableWidget.setItem(i, 2, item)      
                item = QtWidgets.QTableWidgetItem('%.2f'%(self.actual_center_distance))
                self.tableWidget.setItem(i, 3, item)        
                item = QtWidgets.QTableWidgetItem('%.2f'%(speed_ratio_deviation))
                self.tableWidget.setItem(i, 4, item)        
    
                i = i+1  
                
    #=================================================帮助文件========================================================
    @Slot()
    def on_pushButton_second_k1_origin_clicked(self):    
        self.label_help.setPixmap(QtGui.QPixmap("labels/help_k1.png"))
        
    @Slot()
    def on_pushButton_second_k2_origin_clicked(self):    
        self.label_help.setPixmap(QtGui.QPixmap("labels/help_k2.png"))
        
if __name__ == "__main__":
     
    app = QApplication(sys.argv)
    main = MainWindow()
    main.show()
    sys.exit(app.exec_())
