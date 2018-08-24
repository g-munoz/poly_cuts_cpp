Minimize
  0 x(5) + 0 x(6) + 0 x(11) + 0 x(12) - 2 x(13) - 2 x(14) + [ 2 x(13) ^2
   + 2 x(14) ^2 + 2 x(15) ^2 + 2 x(16) ^2 ] / 2 
Subject To
 con1: - 2 x(14) + 2 x(16) - x(8) + x(10) = 0
 con2: - 2 x(13) + 2 x(15) - x(7) + x(9) = 0
 con3: x(16) + x(4) = 1.5
 con4: x(15) + x(3) = 1.5
 con5: - x(16) + x(2) = -0.5
 con6: - x(15) + x(1) = -0.5
 QCP_row_for_con7: [ x(6) * x(12) ] = 0
 QCP_row_for_con8: [ x(7) * x(1) ] = 0
 QCP_row_for_con9: [ x(5) * x(11) ] = 0
 QCP_row_for_con10: [ x(10) * x(4) ] = 0
 QCP_row_for_con11: [ x(8) * x(2) ] = 0
 QCP_row_for_con12: [ x(9) * x(3) ] = 0
Bounds
 x(5) <= 200
 x(6) <= 200
 x(11) <= 200
 x(12) <= 200
 x(8) <= 200
 x(10) <= 200
 x(7) <= 200
 x(9) <= 200
 x(4) <= 200
 x(3) <= 200
 x(2) <= 200
 x(1) <= 200
End
