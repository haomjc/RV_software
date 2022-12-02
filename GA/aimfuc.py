# -*- coding: utf-8 -*-
import numpy as np
from GA.Ui_GA import Ui_Dialog_GA

class aimfunction(Ui_Dialog_GA):

    def __init__(self, parent=None):
        super(aimfunction, self).__init__(parent)
        self.setupUi(self)

    def aimfuc(self, Chrom,LegV):
    
        P = self.doubleSpinBox_P.value()   #输入功率： w
        n = self.doubleSpinBox_n.value()    #转速： rpm 
        n_p = self.spinBox_n_p.value()   
        T_2 = self.doubleSpinBox_T_2.value()    #输出扭矩： N/mm  
        
        #=============================================输入值============================================= 
        Y_1_max = 1.487 
        K_2_min = 1.15
        K_2_max = 2
            
        delta_H_pie = 1100   #摆线轮接触疲劳强度
           
        K = 1.984
        Y_fa = 2.44
        Y_sa = 1.646
        delta_F = 200   #中心轮弯曲疲劳强度
        
        Z_H = 2.5
        Z_E = 189.8
        delta_H = 1100  #渐开线齿轮接触疲劳强度
        
        delta_c = 1100  #滚子接触疲劳强度
        #c_D_o = 0.8  #外滚道直径最大值修正系数
        #摆线轮
        miu = 0.31 #泊松比
        E = 212000  #弹性模量(MPa)
        
        D_o = 25 #行星齿轮轴的直径
    
        #===========================目标变量======================================
        
        b = Chrom[:, [0]]    
        D_z = Chrom[:, [1]]
        d_z = Chrom[:, [2]]
        B = Chrom[:, [3]]
        K_1 = Chrom[:, [4]]
        D_m = Chrom[:, [5]]
        D_r = Chrom[:, [6]]
        L = Chrom[:, [7]]
        
        z_1 = Chrom[:, [8]]
        z_2 = Chrom[:, [9]]
        m_4 = Chrom[:, [10]]
        z_g_2 = Chrom[:, [11]]
        Z = Chrom[:, [12]]
        
        m = m_4/4   #标准模数
        z_g= z_g_2*2+1   #摆线轮齿数为奇数
        #============================================目标函数================================================
        
        T_g = T_2*0.55    
        eta_B1 = 0.99
        eta_B2 = 0.98
        eta_B3 = 0.98    
        eta_B = eta_B1*eta_B2*eta_B3
        
        alpha = np.pi/9
        f = 0.07
        miu_pie = 0.06
        
        F_t = 2*(9550*P/n)/(n_p*m*z_1)   #？
        
        C = 4/(np.pi*z_g)
        i_1 = z_2/z_1
        z_b = z_g+1
        h_a1 = 1*m
        h_a2 = 1*m
        d_a1 = m*z_1+2*h_a1
        d_a2 = m*z_2+2*h_a2
        alpha_a1 = np.arccos(m*z_1*np.cos(alpha)/d_a1)    
        alpha_a2 = np.arccos(m*z_2*np.cos(alpha)/d_a2) 
        epsilon_1 = z_1/(2*np.pi)*(np.tan(alpha_a1)-np.tan(alpha))
        epsilon_2 = z_2/(2*np.pi)*(np.tan(alpha_a2)-np.tan(alpha)) 
        
        i = z_b*z_2/z_1+1
        #==========效率函数=============
        eta_6_62 = 1-f*C/K_1*(1-d_z/D_z)
        eta_1_6 = 1-np.pi*miu_pie*(1/z_1+1/z_2)*(epsilon_1**2+epsilon_2**2+1-epsilon_1-epsilon_2)
        eta_16 = (z_b/z_g-1)*(z_b/z_g-eta_6_62+i_1*z_b/z_g*eta_1_6)/((z_b/z_g-eta_6_62)*(z_b/z_g-1+(i_1*z_b/z_g)))
        eta = -eta_16*eta_B #效率的相反数   
        #=============体积函数===============
        V = np.pi/4*((m*z_1)**2+n_p*(m*z_2)**2)*b+np.pi/4*(D_z+d_z)**2*2*B #体积
        #V = np.pi/4*(D_z+d_z)**2*(2*B+b)
        
        #===============刚度函数====================
        i_16_5 = z_b*z_2/z_1+1
        i_15_6 = z_b*z_2/z_1
        i_65_1 = i_15_6/i_16_5    
        
        #渐开线刚度
        q = 0.04723+0.15551/z_1+0.25791/z_2
        c_th_pie = 1/q
        C_M = 0.8
        C_R = 1
        h_f  = 1.25*m
        C_B = (1+0.5*(1.2-h_f/m))
        epsilon_alpha = epsilon_1+epsilon_2
        c_pie = c_th_pie*C_M*C_R*C_B*b
        c_r = (0.75*epsilon_alpha+0.25)*c_pie*1000  #N/m化为N/mm
        r_s = m*z_1/2
        theta_sp = F_t/(c_r*r_s)
        sigma = theta_sp/i_16_5
        
        #摆线轮刚度    ===================刚度过大====================
        _lambda=0.65
        #计算c_bz
        c_bz = 0
        for phi_i in np.arange(0, np.pi, np.pi/1000):
            S = 1+K_1**2-2*K_1*np.cos(phi_i)
            T = K_1*(2+z_g)*np.cos(phi_i)-(1+(z_g+1)*K_1**2)
            #c_st = np.pi*B*E*D_z*S**(3/2)/(4*(1-miu**2)*(D_z*S**(3/2)+d_z*T+d_z*abs(T)))
            c_st = np.pi*B*E*D_z*S**(3/2)/(4*(1-miu**2)*(D_z*S**(3/2)+2*d_z*abs(T)))            
            e = K_1*D_z/(2*z_b)
            L_i_pie = e*z_g*np.sin(phi_i)/(S**(1/2))
            c_bz = c_bz+_lambda*c_st*L_i_pie**2*(z_b/2)   
            #c_bz = c_bz+_lambda*c_st*(L_i_pie**2)*(z_b/2)
        c_bz = c_bz/1000    #求平均扭转刚度
        beta_H = T_g/c_bz
        beta = i_65_1*z_b/z_g*beta_H    
        
        #轴承刚度
        K_y = 2/np.pi*(1/K_1+(K_1**2-1)/(2*K_1**2)*np.log((1+K_1)/(1-K_1)))
        R_t = T_2*z_2*z_b/(6*i*z_1*K_1*z_b)
        R_x = T_2*z_2*z_b/(3*i*z_1*K_1*D_z)
        R_y = K_y*T_2*z_2*z_b/(3*i*z_1*K_1*D_z)
        R_abs = abs(R_t)+(R_x**2+R_y**2)**(1/2)
        
        K_r = 0.340*10**4*(R_abs**0.1*Z**0.9*L**0.8)   
        
        delta_cb = R_abs/K_r
        a_sp = m*(z_1+z_2)/2
        tau = z_g*delta_cb*i_65_1/(z_b*a_sp)
        
        theta = sigma+tau+beta
        K_pie = -T_2/theta  #刚度的相反数(N·mm/rad)
        
        # ==================================约束条件============================================
    
        #行星轮齿比关系：
        #idx1 = np.where(1.5*z_1 > z_2)[0]
        #不根切条件：
        #idx2 = np.where(z_1  < 17)[0]
        #最小模数：
        #idx3 = np.where(m < 1)[0]
        #摆线齿廓不根切:
        idx4 = np.where((d_z/D_z) >((27*z_g*(1-K_1**2))**(1/2)/((z_g+2)**(3/2))))[0]    
        # 针齿分布圆直径:
        #idx5 = np.where(D_z<21*T_2**(1/3))[0]
        #idx6 = np.where(D_z>26*T_2**(1/3))[0]    
        #摆线轮接触强度:
        idx7 = np.where(759*(T_g*Y_1_max/(B*D_z**2))**(1/2)>delta_H_pie)[0]  #注意T_g 单位换算
        #针径系数:
    
        K_2 = D_z*np.sin(np.pi/z_g)/d_z
        idx8 = np.where(K_2<K_2_min)[0]
        idx9 = np.where(K_2>K_2_max)[0]    
        #摆线轮宽度：
        #idx10 = np.where(B<0.05*D_z)[0]
        #idx11 = np.where(B>0.1*D_z)[0]    
        #行星轮模数与齿宽的关系： 
        #idx12 = np.where(b<5*m)[0]
        #idx13 = np.where(b>17*m)[0]      
        
        #一级传动与二级传动大小一致：
        #idx14 = np.where(m*z_1+2*m*z_2<0.8*(D_z+d_z))[0]
        idx15 = np.where(m*z_1+2*m*z_2>1.1*(D_z+d_z))[0] 
        #邻接条件:
        idx16 = np.where((z_1+z_2)*np.sin(np.pi/n_p)<z_2+2)[0]
        #中心轮弯曲疲劳强度：
        idx17 = np.where(K*F_t*Y_fa*Y_sa/(b*m)>delta_F)[0]
        #渐开线齿轮接触疲劳强度：
        u = z_2/z_1
        idx18 = np.where((K*F_t/(b*m*z_1)*(u+1)/u)**(1/2)*Z_H*Z_E>delta_H)[0]
        #轴承寿命： 
        n_b = n/i*z_b
        p = R_abs*(((R_x**2+R_y**2)**(1/2)/R_abs-0.5)**2+0.75)   
        b_m = 1.1
        lambda_nu = 0.83
        gamma = D_r/D_m
        f_c = 208*lambda_nu*gamma**(2/9)*(1-gamma)**(29/27)/((1+gamma)**(1/4))*(1+(1.04*((1-gamma)/(1+gamma))**(143/108))**(9/2))**(-2/9)
        C_d = b_m*f_c*L**(7/9)*Z**(3/4)*D_r**(29/27)
        L_10h = 10**6/(60*n_b)*(C_d/p)**(10/3)
        idx19 = np.where(L_10h<10000)[0]   #？
        #圆柱滚子轴承的接触长度：
        idx20 = np.where(L<D_r)[0]
        idx21 = np.where(L>2.5*D_r)[0]
        #滚子直径：
        Q_max = 4.08/1000*R_abs/Z
        #D_o = m*z_2*c_D_o  #外滚道直径最大值
        s = 5  #安全系数
        D_r_min = 268.71*Q_max**(1/2)/(delta_c/s)
        
        idx22 = np.where(D_r<D_r_min)[0]
        #idx23 = np.where(D_r>(D_r+D_m)*np.sin(np.pi/5)/(1+np.sin(np.pi/5)))[0]

        idx23 = np.where(D_m-D_r>D_o+2*e)[0]
        #滚子数:
        #idx24 = np.where(Z<5)[0]
        #idx25 = np.where(Z>np.pi/np.arcsin(D_r/D_m))[0]
        #传动比范围:
    
        i_min = T_2/(9550*P/n)/(-eta)
        i_max = i_min*1.2
        
        idx26 = np.where(i<i_min)[0]
        idx27 = np.where(i>i_max)[0]
        #轴承节圆直径:
        idx28 = np.where(D_m>m*z_2-D_r)[0]
        
        #圆柱滚子轴承的接触长度：
        idx29 = np.where(L>B)[0]
        
        #行星齿轮装配条件
        #if (n_p  == 2 and (z_2 - z_1)%2 ==0) or (n_p  == 3 and z_2%2 == 0 and (2*z_1 - z_2)%6 == 0) or (n_p  == 3 and z_2%2 == 1 and (2*z_1 - z_2 )%6 == 3):
        if n_p == 2:
            asb1 = (z_2 - z_1)%2
            idx30 = np.where(z_1*asb1 != 0)[0]
        
        elif n_p ==3:
            asb2 = z_2%2
            asb3 = (2*z_1 - z_2)%6
            idx30 = np.where(z_1*(asb2*3-asb3) != 0)[0]
            
        #避免滚子相碰：
        idx31 = np.where(2*Z*np.arcsin(D_r/D_m)+Z*(np.pi/180)>2*np.pi)[0]
        
        #=======================目标函数的限制条件===============================
        #idx32 = np.where(-eta<0.84)[0]
        #idx33 = np.where(V>1800000)[0]

        exIdx = np.unique(np.hstack([idx4,idx7,idx8,idx9, idx15, idx16,idx17,idx18, idx19,idx20, idx21,idx22,idx23, idx26, idx27, idx28, idx29, idx30, idx31])) # 得到非可行解的下标
        LegV[exIdx] = 0 # 标记非可行解对应的可行性列向量中元素的值为0
        return [np.hstack([eta, V, K_pie]), LegV]

