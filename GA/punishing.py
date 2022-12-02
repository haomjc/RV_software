# -*- coding: utf-8 -*-
import numpy as np

def punishing(LegV, FitnV):
    FitnV[np.where(LegV == 0)[0]] = np.min(FitnV)//2  # 对非可行解惩罚
    return FitnV
