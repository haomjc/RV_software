# -*- coding: utf-8 -*-
import numpy as np
import geatpy as ga # 导入geatpy库
import time
from GA.Ui_GA import Ui_Dialog_GA
from PySide2.QtWidgets import QApplication

class nsga2(Ui_Dialog_GA):
    
    def moea_nsga2(self, AIM_M, AIM_F, PUN_M, PUN_F, FieldDR, problem,R_num,  maxormin, MAXGEN, MAXSIZE, NIND, SUBPOP, GGAP, selectStyle, recombinStyle, recopt, pm, distribute):
        
        # 获取目标函数和罚函数
        aimfuc = getattr(AIM_M, AIM_F) # 获得目标函数
        if PUN_F is not None:
            punishing = getattr(PUN_M, PUN_F) # 获得罚函数
        #==========================初始化配置===========================
    
        # 获取目标函数和罚函数
        aimfuc = getattr(AIM_M, AIM_F) # 获得目标函数
        #=========================开始遗传算法进化=======================
        if problem == 'R':
            Chrom = ga.crtrp(NIND, FieldDR) # 生成实数值种群
        elif problem == 'I':
            Chrom = ga.crtip(NIND, FieldDR) # 生成整数值种群
        elif problem == 'M':      #生成混合种群
            Chrom = np.hstack([ga.crtrp(NIND, FieldDR[:, 0:R_num]), ga.crtip(NIND, FieldDR[:,R_num: ])])   
    
        LegV = np.ones((NIND, 1)) # 初始化可行性列向量
        [ObjV, LegV] = aimfuc.aimfunction.aimfuc(self, Chrom, LegV) # 计算种群目标函数值
        NDSet = np.zeros((0, Chrom.shape[1])) # 定义帕累托最优解集合(初始为空集)
        NDSetObjV = np.zeros((0, ObjV.shape[1])) # 定义帕累托最优解对应的目标函数集合(初始为空集)
        ax = None # 存储上一桢动画
        start_time = time.time() # 开始计时
        # 计算初代
        [FitnV, levels] = ga.ndomindeb(maxormin * ObjV, 1, LegV) # deb非支配分级
        frontIdx = np.where(levels == 1)[0] # 处在第一级的个体即为种群的非支配个体
        if PUN_F is not None:
            FitnV = punishing.punishing(LegV, FitnV) # 调用罚函数
        # 更新帕累托最优集以及种群非支配个体的适应度
        [FitnV, NDSet, NDSetObjV, repnum] = ga.upNDSet(Chrom, maxormin * ObjV, FitnV, NDSet, maxormin * NDSetObjV, frontIdx, LegV)
        NDSetObjV *= maxormin # 还原在传入upNDSet函数前被最小化处理过的NDSetObjV
        [NDSet, NDSetObjV] = ga.redisNDSet(NDSet, NDSetObjV, NDSetObjV.shape[1] * MAXSIZE) # 利用拥挤距离选择帕累托前沿的子集，在进化过程中最好比上限多筛选出几倍的点集
        # 开始进化！！
        for gen in range(MAXGEN):
            # 进行遗传操作！！
            SelCh=ga.recombin(recombinStyle, Chrom, recopt, SUBPOP) #交叉
            if problem == 'R':
                SelCh=ga.mutbga(SelCh,FieldDR, pm) # 变异
                if repnum >= Chrom.shape[0] * 0.01: # 当最优个体重复率高达1%时，进行一次高斯变异
                    SelCh=ga.mutgau(SelCh, FieldDR, pm) # 高斯变异
            elif problem == 'I':
                SelCh=ga.mutint(SelCh, FieldDR, pm)
            elif problem == 'M':
                SelCh_R=ga.mutbga(SelCh[:, 0:R_num],FieldDR[:, 0:R_num], pm) # 变异
                if repnum >= Chrom.shape[0] * 0.01: # 当最优个体重复率高达1%时，进行一次高斯变异
                    SelCh_R=ga.mutgau(SelCh[:, 0:R_num], FieldDR[:, 0:R_num], pm) # 高斯变异
                SelCh_I=ga.mutint(SelCh[:,R_num: ], FieldDR[:,R_num: ], pm)     
                SelCh= np.hstack([SelCh_R,SelCh_I]) 
                    
            [ObjVSel, LegVSel] = aimfuc.aimfunction.aimfuc(self, SelCh, LegV) # 求育种个体的目标函数值
            # 父子合并
            Chrom = np.vstack([Chrom, SelCh])
            ObjV = np.vstack([ObjV, ObjVSel])
            LegV = np.vstack([LegV, LegVSel])
            [FitnV, levels] = ga.ndomindeb(maxormin * ObjV, 1, LegV) # deb非支配分级
            frontIdx = np.where(levels == 1)[0] # 处在第一级的个体即为种群的非支配个体
            if PUN_F is not None:
                FitnV = punishing.punishing(LegV, FitnV) # 调用罚函数
            # 更新帕累托最优集以及种群非支配个体的适应度
            [FitnV, NDSet, NDSetObjV, repnum] = ga.upNDSet(Chrom, maxormin * ObjV, FitnV, NDSet, maxormin * NDSetObjV, frontIdx, LegV)
            NDSetObjV *= maxormin # 还原在传入upNDSet函数前被最小化处理过的NDSetObjV
            [NDSet, NDSetObjV] = ga.redisNDSet(NDSet, NDSetObjV, NDSetObjV.shape[1] * MAXSIZE) # 利用拥挤距离选择帕累托前沿的子集，在进化过程中最好比上限多筛选出几倍的点集
            if distribute == True: # 若要增强种群的分布性(可能会导致帕累托前沿搜索效率降低)
                # 计算每个目标下相邻个体的距离(不需要严格计算欧氏距离)
                for i in range(ObjV.shape[1]):
                    idx = np.argsort(ObjV[:, i], 0)
                    dis = np.diff(ObjV[idx, i]) / (np.max(ObjV[idx, i]) - np.min(ObjV[idx, i]) + 1) # 差分计算距离的偏移量占比，即偏移量除以目标函数的极差。加1是为了避免极差为0
                    dis = np.hstack([dis, dis[-1]])
                    FitnV[idx, 0] *= np.exp(dis) # 根据相邻距离修改适应度，突出相邻距离大的个体，以增加种群的多样性
            [Chrom, ObjV, LegV] = ga.selecting(selectStyle, Chrom, FitnV, GGAP, SUBPOP, ObjV, LegV) # 选择出下一代
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
        print(NDSet)
        print('单位时间找到帕累托前沿点个数：%s 个'%(int(NDSet.shape[0] // times)))
        '''
        # 返回帕累托最优集以及执行时间
        return [ObjV, NDSet, NDSetObjV, end_time - start_time]
