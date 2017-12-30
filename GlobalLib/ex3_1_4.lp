Minimize
  - 2 x(1) + x(2) - x(3)
Subject To
 con1: 3 x(2) + x(3) <= 6
 con2: x(1) + x(2) + x(3) <= 4
 QCP_row_for_con3: - 20 x(1) + 9 x(2) - 13 x(3) + [ 4 x(1) ^ 2 - 4 x(1) * x(2)
   + 4 x(1) * x(3) + 2 x(2) ^ 2 - 2 x(2) * x(3) + 2 x(3) ^ 2 ] >= -24
Bounds
 x(1) <= 2
 x(2) <= 2
 x(3) <= 3
End
