import numpy as np

P = 860   #输入功率： w
n = 1880    #转速： rpm 
n_p = 3
T_2 = 412000    #输出扭矩： N/mm  

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

b = 10
D_z = 165
d_z = 7
B = 14
K_1 = 0.7273
D_m = 41.2
D_r = 8
L = 12

z_1 = 19
z_2 = 57


Z = 16

m = 1.25   #标准模数
z_g= 39   #摆线轮齿数为奇数
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

print(eta, V, K_pie, i_16_5)
