% 函数功能：检测两条线段是否相交，后续问题需要判断是否相交

function is_intersecting = check_intersection(P1, P2, Q1, Q2)  

    cross_product = @(A, B) A(1)*B(2) - A(2)*B(1);
    
    P1P2 = P2 - P1;
    P1Q1 = Q1 - P1;
    P1Q2 = Q2 - P1;
    
    Q1Q2 = Q2 - Q1;
    Q1P1 = P1 - Q1;
    Q1P2 = P2 - Q1;
    
    d1 = cross_product(P1P2, P1Q1);
    d2 = cross_product(P1P2, P1Q2);
    d3 = cross_product(Q1Q2, Q1P1);
    d4 = cross_product(Q1Q2, Q1P2);
    
    % 如果符号不同，表示相交
    is_intersecting = (d1 * d2 < 0) && (d3 * d4 < 0);  
end
