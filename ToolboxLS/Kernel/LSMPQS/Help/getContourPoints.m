 function [C1 total] = getContourPoints(C)
 
 % function [C1 total] = getContourPoints(C)
 % C = matrix output from contourc or contour functions, might have > 1
 % contiguous drawing segments which include a column with value of 
 % the contour and number of points and contour points (x,y) in columns
 % C1 = has only (x,y) pairs w/o any other information
 % total - number of columns in C1
 
 [tmp N ] = size(C);  %note tmp always ==2
    
 total = 0;
 i = 1;
 n = C(2,1);
 total = total + n;
 s = 2; e = n+1;
 C1 = C(:,s:e); 
 while (e < N)
     i = i+1;
     n = C(2,e+1);
     total = total + n;
     s = e + 2;
     e = e + n + 1;
     [tmp e0] = size(C1);
     C1(:,e0+1:e0+n) = C(:,s:e);
 end