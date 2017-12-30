Minimize
  0 x(1) - 9 x(4) - 15 x(5) + 6 x(6) + 16 x(7) + 10 x(8) + 10 x(9)
Subject To
 con1: x(5) - x(9) - x(3) = 0
 con2: x(4) - x(8) - x(2) = 0
 con3: - x(6) - x(7) + x(3) + x(2) = 0
 QCP_row_for_con4: - 3 x(6) - x(7) + [ x(1) * x(3) + x(1) * x(2) ] = 0
 QCP_row_for_con5: - 1.5 x(5) + 2 x(9) + [ x(1) * x(3) ] <= 0
 QCP_row_for_con6: - 2.5 x(4) + 2 x(8) + [ x(1) * x(2) ] <= 0
Bounds
 x(1) <= 500
 x(4) <= 100
 x(5) <= 200
 x(6) <= 500
 x(7) <= 500
 x(8) <= 500
 x(9) <= 500
 x(3) <= 500
 x(2) <= 500
End
