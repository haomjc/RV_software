# -*- coding: utf-8 -*-

import math
from PySide2.QtCore import Slot
from PySide2.QtWidgets import QDialog
from PySide2 import QtWidgets

from Ui_teethmatch_tool_window import Ui_Dialog_teethmatch_tool


class Dialog_teethmatch_tool(QDialog, Ui_Dialog_teethmatch_tool):

    def __init__(self, parent=None):

        super(Dialog_teethmatch_tool, self).__init__(parent)
        self.setupUi(self)
        self.radioButton_actual_center_distance = False
        self.radioButton_normal_module = False
        self.speed_ratio_deviation = 0.05    
        
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
    def on_pushButton_calculate_clicked(self):

        self.speed_ratio = self.doubleSpinBox_speed_ratio.value()
        self.gearteeth_min = self.spinBox_gearteeth_min.value()
        self.beta_angle = self.doubleSpinBox_beta_angle.value()/180*math.pi
        self.d_or_mn = self.doubleSpinBox_d_or_mn.value()
        
        if self.radioButton_actual_center_distance == True:
            self.actual_center_distance = self.d_or_mn
        elif self.radioButton_normal_module == True:
            self.normal_module = self.d_or_mn
            
        i = 0
        self.tableWidget.setRowCount(int((1+self.speed_ratio_deviation)*self.gearteeth_min*self.speed_ratio)-math.ceil((1-self.speed_ratio_deviation)*self.gearteeth_min*self.speed_ratio))

        for gearteeth_max in range(math.ceil((1-self.speed_ratio_deviation)*self.gearteeth_min*self.speed_ratio),int((1+self.speed_ratio_deviation)*self.gearteeth_min*self.speed_ratio) ):

            speed_ratio_deviation = gearteeth_max/self.gearteeth_min
            if self.radioButton_actual_center_distance == True:
                self.normal_module = 2*self.actual_center_distance*math.cos(self.beta_angle)/(self.gearteeth_min*(speed_ratio_deviation+1))
            elif self.radioButton_normal_module == True:
                self.actual_center_distance = self.gearteeth_min*self.normal_module*(speed_ratio_deviation+1)/(2*math.cos(self.beta_angle))

            item = QtWidgets.QTableWidgetItem('%d'% (self.gearteeth_min))
            self.tableWidget.setItem(i, 0, item)            
            item = QtWidgets.QTableWidgetItem('%d'%(gearteeth_max))
            self.tableWidget.setItem(i, 1, item)                        
            item = QtWidgets.QTableWidgetItem('%.2f'%(self.normal_module))
            self.tableWidget.setItem(i, 2, item) 
            item = QtWidgets.QTableWidgetItem('%.2f'%(self.beta_angle/math.pi*180))
            self.tableWidget.setItem(i, 3, item)       
            item = QtWidgets.QTableWidgetItem('%.2f'%(self.actual_center_distance))
            self.tableWidget.setItem(i, 4, item)        
            item = QtWidgets.QTableWidgetItem('%.2f'%(speed_ratio_deviation))
            self.tableWidget.setItem(i, 5, item)        

            i = i+1
            
    @Slot()
    def on_pushButton_ok_clicked(self):
        Dialog_teethmatch_tool.hide(self)

