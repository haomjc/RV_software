# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'F:\OneDrive\PythonFile\RV_software\RV_1203\K_a_coefficient.ui',
# licensing of 'F:\OneDrive\PythonFile\RV_software\RV_1203\K_a_coefficient.ui' applies.
#
# Created: Tue Dec  4 14:11:08 2018
#      by: pyside2-uic  running on PySide2 5.11.2
#
# WARNING! All changes made in this file will be lost!

from PySide2 import QtCore, QtGui, QtWidgets

class Ui_Dialog(object):
    def setupUi(self, Dialog):
        Dialog.setObjectName("Dialog")
        Dialog.resize(555, 297)
        Dialog.setSizeGripEnabled(True)
        self.tableWidget = QtWidgets.QTableWidget(Dialog)
        self.tableWidget.setGeometry(QtCore.QRect(60, 60, 461, 151))
        self.tableWidget.setObjectName("tableWidget")
        self.tableWidget.setColumnCount(4)
        self.tableWidget.setRowCount(4)
        item = QtWidgets.QTableWidgetItem()
        self.tableWidget.setVerticalHeaderItem(0, item)
        item = QtWidgets.QTableWidgetItem()
        self.tableWidget.setVerticalHeaderItem(1, item)
        item = QtWidgets.QTableWidgetItem()
        self.tableWidget.setVerticalHeaderItem(2, item)
        item = QtWidgets.QTableWidgetItem()
        self.tableWidget.setVerticalHeaderItem(3, item)
        item = QtWidgets.QTableWidgetItem()
        self.tableWidget.setHorizontalHeaderItem(0, item)
        item = QtWidgets.QTableWidgetItem()
        self.tableWidget.setHorizontalHeaderItem(1, item)
        item = QtWidgets.QTableWidgetItem()
        self.tableWidget.setHorizontalHeaderItem(2, item)
        item = QtWidgets.QTableWidgetItem()
        self.tableWidget.setHorizontalHeaderItem(3, item)
        item = QtWidgets.QTableWidgetItem()
        self.tableWidget.setItem(0, 0, item)
        item = QtWidgets.QTableWidgetItem()
        self.tableWidget.setItem(0, 1, item)
        item = QtWidgets.QTableWidgetItem()
        self.tableWidget.setItem(0, 2, item)
        item = QtWidgets.QTableWidgetItem()
        self.tableWidget.setItem(0, 3, item)
        item = QtWidgets.QTableWidgetItem()
        self.tableWidget.setItem(1, 0, item)
        item = QtWidgets.QTableWidgetItem()
        self.tableWidget.setItem(1, 1, item)
        item = QtWidgets.QTableWidgetItem()
        self.tableWidget.setItem(1, 2, item)
        item = QtWidgets.QTableWidgetItem()
        self.tableWidget.setItem(1, 3, item)
        item = QtWidgets.QTableWidgetItem()
        self.tableWidget.setItem(2, 0, item)
        item = QtWidgets.QTableWidgetItem()
        self.tableWidget.setItem(2, 1, item)
        item = QtWidgets.QTableWidgetItem()
        self.tableWidget.setItem(2, 2, item)
        item = QtWidgets.QTableWidgetItem()
        self.tableWidget.setItem(2, 3, item)
        item = QtWidgets.QTableWidgetItem()
        self.tableWidget.setItem(3, 0, item)
        item = QtWidgets.QTableWidgetItem()
        self.tableWidget.setItem(3, 1, item)
        item = QtWidgets.QTableWidgetItem()
        self.tableWidget.setItem(3, 2, item)
        item = QtWidgets.QTableWidgetItem()
        self.tableWidget.setItem(3, 3, item)
        self.label = QtWidgets.QLabel(Dialog)
        self.label.setGeometry(QtCore.QRect(10, 10, 111, 31))
        self.label.setObjectName("label")
        self.verticalLayoutWidget = QtWidgets.QWidget(Dialog)
        self.verticalLayoutWidget.setGeometry(QtCore.QRect(70, 220, 386, 61))
        self.verticalLayoutWidget.setObjectName("verticalLayoutWidget")
        self.verticalLayout = QtWidgets.QVBoxLayout(self.verticalLayoutWidget)
        self.verticalLayout.setContentsMargins(0, 0, 0, 0)
        self.verticalLayout.setObjectName("verticalLayout")
        self.label_2 = QtWidgets.QLabel(self.verticalLayoutWidget)
        self.label_2.setObjectName("label_2")
        self.verticalLayout.addWidget(self.label_2)
        self.label_3 = QtWidgets.QLabel(self.verticalLayoutWidget)
        self.label_3.setObjectName("label_3")
        self.verticalLayout.addWidget(self.label_3)
        self.label_4 = QtWidgets.QLabel(self.verticalLayoutWidget)
        self.label_4.setObjectName("label_4")
        self.verticalLayout.addWidget(self.label_4)
        self.label_5 = QtWidgets.QLabel(Dialog)
        self.label_5.setGeometry(QtCore.QRect(10, 120, 41, 20))
        self.label_5.setObjectName("label_5")
        self.label_6 = QtWidgets.QLabel(Dialog)
        self.label_6.setGeometry(QtCore.QRect(270, 30, 91, 16))
        self.label_6.setObjectName("label_6")
        self.label_7 = QtWidgets.QLabel(Dialog)
        self.label_7.setGeometry(QtCore.QRect(0, 140, 54, 20))
        self.label_7.setObjectName("label_7")

        self.retranslateUi(Dialog)
        QtCore.QMetaObject.connectSlotsByName(Dialog)

    def retranslateUi(self, Dialog):
        Dialog.setWindowTitle(QtWidgets.QApplication.translate("Dialog", "Dialog", None, -1))
        self.tableWidget.verticalHeaderItem(0).setText(QtWidgets.QApplication.translate("Dialog", "均匀平稳", None, -1))
        self.tableWidget.verticalHeaderItem(1).setText(QtWidgets.QApplication.translate("Dialog", "轻微冲击", None, -1))
        self.tableWidget.verticalHeaderItem(2).setText(QtWidgets.QApplication.translate("Dialog", "中等冲击", None, -1))
        self.tableWidget.verticalHeaderItem(3).setText(QtWidgets.QApplication.translate("Dialog", "严重冲击", None, -1))
        self.tableWidget.horizontalHeaderItem(0).setText(QtWidgets.QApplication.translate("Dialog", "均匀平稳", None, -1))
        self.tableWidget.horizontalHeaderItem(1).setText(QtWidgets.QApplication.translate("Dialog", "轻微冲击", None, -1))
        self.tableWidget.horizontalHeaderItem(2).setText(QtWidgets.QApplication.translate("Dialog", "中等冲击", None, -1))
        self.tableWidget.horizontalHeaderItem(3).setText(QtWidgets.QApplication.translate("Dialog", "严重冲击", None, -1))
        __sortingEnabled = self.tableWidget.isSortingEnabled()
        self.tableWidget.setSortingEnabled(False)
        self.tableWidget.item(0, 0).setText(QtWidgets.QApplication.translate("Dialog", "1.00\n"
"", None, -1))
        self.tableWidget.item(0, 1).setText(QtWidgets.QApplication.translate("Dialog", "1.25", None, -1))
        self.tableWidget.item(0, 2).setText(QtWidgets.QApplication.translate("Dialog", "1.50", None, -1))
        self.tableWidget.item(0, 3).setText(QtWidgets.QApplication.translate("Dialog", "1.75", None, -1))
        self.tableWidget.item(1, 0).setText(QtWidgets.QApplication.translate("Dialog", "1.10", None, -1))
        self.tableWidget.item(1, 1).setText(QtWidgets.QApplication.translate("Dialog", "1.35", None, -1))
        self.tableWidget.item(1, 2).setText(QtWidgets.QApplication.translate("Dialog", "1.60", None, -1))
        self.tableWidget.item(1, 3).setText(QtWidgets.QApplication.translate("Dialog", "1.85", None, -1))
        self.tableWidget.item(2, 0).setText(QtWidgets.QApplication.translate("Dialog", "1.25", None, -1))
        self.tableWidget.item(2, 1).setText(QtWidgets.QApplication.translate("Dialog", "1.50", None, -1))
        self.tableWidget.item(2, 2).setText(QtWidgets.QApplication.translate("Dialog", "1.75", None, -1))
        self.tableWidget.item(2, 3).setText(QtWidgets.QApplication.translate("Dialog", "2.00", None, -1))
        self.tableWidget.item(3, 0).setText(QtWidgets.QApplication.translate("Dialog", "1.50", None, -1))
        self.tableWidget.item(3, 1).setText(QtWidgets.QApplication.translate("Dialog", "1.75", None, -1))
        self.tableWidget.item(3, 2).setText(QtWidgets.QApplication.translate("Dialog", "2.00", None, -1))
        self.tableWidget.item(3, 3).setText(QtWidgets.QApplication.translate("Dialog", "2.25或更大", None, -1))
        self.tableWidget.setSortingEnabled(__sortingEnabled)
        self.label.setText(QtWidgets.QApplication.translate("Dialog", "使用系数KA", None, -1))
        self.label_2.setText(QtWidgets.QApplication.translate("Dialog", "注：1.对于增速传动，根据经验建议取表中值的1.1倍。", None, -1))
        self.label_3.setText(QtWidgets.QApplication.translate("Dialog", "    2.当外部机械与齿轮装置之间有挠性联接时，通常KA值可适当减小。", None, -1))
        self.label_4.setText(QtWidgets.QApplication.translate("Dialog", "    3.在采用表荐值时，至少应取最小弯曲强度安全系数SFmin=1.25。", None, -1))
        self.label_5.setText(QtWidgets.QApplication.translate("Dialog", "原动机", None, -1))
        self.label_6.setText(QtWidgets.QApplication.translate("Dialog", "工作机工作特性", None, -1))
        self.label_7.setText(QtWidgets.QApplication.translate("Dialog", "工作特性", None, -1))


if __name__ == "__main__":
    import sys
    app = QtWidgets.QApplication(sys.argv)
    Dialog = QtWidgets.QDialog()
    ui = Ui_Dialog()
    ui.setupUi(Dialog)
    Dialog.show()
    sys.exit(app.exec_())
