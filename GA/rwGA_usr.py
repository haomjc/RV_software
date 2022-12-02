# -*- coding: utf-8 -*-
import numpy as np
import geatpy as ga # 导入geatpy库
import time
from GA.Ui_GA import Ui_Dialog_GA
from PySide2.QtWidgets import QApplication

class rwGA(Ui_Dialog_GA):

    def moea_rwGA(self,AIM_M, AIM_F, PUN_M, PUN_F, FieldDR, problem, R_num,maxormin, MAXGEN, MAXSIZE, NIND, SUBPOP, GGAP, selectStyle, recombinStyle, recopt, pm, distribute, drawing = 1):
        
        #==========================初始化配置===========================
        # 获取目标函数和罚函数
        aimfuc = getattr(AIM_M, AIM_F) # 获得目标函数
        if PUN_F is not None:
            punishing = getattr(PUN_M, PUN_F) # 获得罚函数
    
        #=========================开始遗传算法进化=======================
        if problem == 'R':
            Chrom = ga.crtrp(NIND, FieldDR) # 生成实数值种群
        elif problem == 'I':
            Chrom = ga.crtip(NIND, FieldDR) # 生成整数值种群
        elif problem == 'M':      #生成混合种群
            Chrom = np.hstack([ga.crtrp(NIND, FieldDR[:, 0:R_num]), ga.crtip(NIND, FieldDR[:,R_num: ])])           
            
        LegV = np.ones((NIND, 1)) # 初始化种群的可行性列向量
        [ObjV, LegV] = aimfuc.aimfunction.aimfuc(self,Chrom, LegV) # 计算种群目标函数值
        NDSet = np.zeros((0, Chrom.shape[1])) # 定义帕累托最优解记录器
        NDSetObjV = np.zeros((0, ObjV.shape[1])) # 定义帕累托最优解的目标函数值记录器
        ax = None # 存储上一桢动画
        start_time = time.time() # 开始计时
        # 开始进化！！
        for gen in range(MAXGEN):
            [CombinObjV, weight] = ga.rwGA(maxormin * ObjV, LegV) # 计算适应性权重以及多目标的加权单目标
            CombinObjV *= maxormin # 还原在传入函数前被最小化处理过的目标函数值
            FitnV  = ga.ranking(maxormin * CombinObjV, LegV, None, SUBPOP) # 根据加权单目标计算适应度
            if PUN_F is not None:
                FitnV = punishing.punishing(LegV, FitnV) # 调用罚函数作进一步的惩罚
            # 更新帕累托最优集以及种群非支配个体的适应度
            [FitnV, NDSet, NDSetObjV, repnum] = ga.upNDSet(Chrom, maxormin * ObjV, FitnV, NDSet, maxormin * NDSetObjV, None, LegV)
            NDSetObjV *= maxormin # 还原在传入upNDSet函数前被最小化处理过的NDSetObjV
            [NDSet, NDSetObjV] = ga.redisNDSet(NDSet, NDSetObjV, NDSetObjV.shape[1] * MAXSIZE) # 利用拥挤距离选择帕累托前沿的子集，在进化过程中最好比上限多筛选出几倍的点集
            if distribute == True: # 若要增强种群的分布性(可能会导致帕累托前沿搜索效率降低)
                # 计算每个目标下相邻个体的距离(不需要严格计算欧氏距离)
                for i in range(ObjV.shape[1]):
                    idx = np.argsort(ObjV[:, i], 0)
                    dis = np.diff(ObjV[idx, i]) / (np.max(ObjV[idx, i]) - np.min(ObjV[idx, i]) + 1) # 差分计算距离的偏移量占比，即偏移量除以目标函数的极差。加1是为了避免极差为0
                    dis = np.hstack([dis, dis[-1]])
                    FitnV[idx, 0] *= np.exp(dis) # 根据相邻距离修改适应度，突出相邻距离大的个体，以增加种群的多样性
            # 进行遗传操作！！
            SelCh=ga.selecting(selectStyle, Chrom, FitnV, GGAP, SUBPOP) # 选择
            SelCh=ga.recombin(recombinStyle, SelCh, recopt, SUBPOP) #交叉
            if problem == 'R':
                SelCh=ga.mutbga(SelCh,FieldDR, pm) # 变异
                if repnum > Chrom.shape[0] * 0.01: # 当最优个体重复率高达1%时，进行一次高斯变异
                    SelCh=ga.mutgau(SelCh, FieldDR, pm) # 高斯变异
            elif problem == 'I':
                SelCh=ga.mutint(SelCh, FieldDR, pm)
            elif problem == 'M':
                SelCh_R=ga.mutbga(SelCh[:, 0:R_num],FieldDR[:, 0:R_num], pm) # 变异
                if repnum >= Chrom.shape[0] * 0.01: # 当最优个体重复率高达1%时，进行一次高斯变异
                    SelCh_R=ga.mutgau(SelCh[:, 0:R_num], FieldDR[:, 0:R_num], pm) # 高斯变异
                SelCh_I=ga.mutint(SelCh[:,R_num: ], FieldDR[:,R_num: ], pm)     
                SelCh= np.hstack([SelCh_R,SelCh_I])           
                
            LegVSel = np.ones((SelCh.shape[0], 1)) # 初始化育种种群的可行性列向量
            [ObjVSel, LegVSel] = aimfuc.aimfunction.aimfuc(self,SelCh,LegVSel) # 求育种个体的目标函数值
            [CombinObjV, weight] = ga.rwGA(maxormin * ObjVSel, LegVSel)
            CombinObjV *= maxormin # 还原在传入函数前被最小化处理过的目标函数值
            FitnVSel = ga.ranking(maxormin * CombinObjV, LegVSel, None, SUBPOP)
            if PUN_F is not None:
                FitnVSel = punishing.punishing(LegVSel, FitnVSel) # 调用罚函数
            [Chrom,ObjV,LegV] = ga.reins(Chrom,SelCh,SUBPOP,1,0.9,FitnV,FitnVSel,ObjV,ObjVSel,LegV,LegVSel) #重插入
            
            QApplication.processEvents()
            self.progressBar.setProperty("value", (gen+1)/MAXGEN*100)
            
        end_time = time.time() # 结束计时
        [NDSet, NDSetObjV] = ga.redisNDSet(NDSet, NDSetObjV, MAXSIZE) # 最后根据拥挤距离选择均匀分布的点
        #=========================绘图及输出结果=========================
        '''
        if drawing != 0:
            ga.frontplot(NDSetObjV,True)
        times = end_time - start_time
        print('用时：%s 秒'%(times))
        print('帕累托前沿点个数：%s 个'%(NDSet.shape[0]))
        print('单位时间找到帕累托前沿点个数：%s 个'%(int(NDSet.shape[0] // times)))
        '''
        # 返回帕累托最优集以及执行时间
        return [ObjV, NDSet, NDSetObjV, end_time - start_time]
