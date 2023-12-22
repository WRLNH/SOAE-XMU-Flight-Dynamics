# Intro

这是厦门大学航空航天学院飞行器设计与工程专业课程, 飞行力学的代码备份.

# Notes

`UAV_CD.m`, `UAV_CL.m`, `UAV_CM.m`和`UAV_CY.m`这些气动系数函数有别于`uavL.m`内的相应名称函数.

单独出来的气动系数函数传入的角度单位为角度制. `uavL.m`内部的气动系数函数传入的角度单位为弧度制.

`UAV_CD.m`, `UAV_CL.m`, `UAV_CM.m`和`UAV_CY.m`在`My_Matrix.m`有使用, 因此如果对这些文件做出改动，则需要检查是否会对`My_Matrix.m`有影响.

`Flight_Envelope_Analysis.m`尚未完成.
