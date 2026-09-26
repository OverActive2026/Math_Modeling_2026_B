function test_q4_online_rule_controller()
%TEST_Q4_ONLINE_RULE_CONTROLLER Check each rule without slow stack solves.
    healthyV = 0.8*ones(1,5);
    normalRise = 0.6*ones(1,5);
    zeroVrate = zeros(1,5);

    [q,mode] = q4_online_rule_controller( ...
        -30*ones(1,5),healthyV,zeros(1,5),zeroVrate);
    assert(isequal(q,ones(1,5)));
    assert(mode.insufficientRise);

    [q,mode] = q4_online_rule_controller( ...
        [-1.4 13 17 13 -1.4],healthyV,normalRise,zeroVrate);
    assert(q(3) > 0 && q(3) < 1);
    assert(q(2) == 1 && q(4) == 1);
    assert(mode.reducedPower(3));

    [q,~] = q4_online_rule_controller( ...
        [-1.1 13 16 13 -1.1],healthyV,normalRise,zeroVrate);
    assert(q(3) == 0);

    [q,mode] = q4_online_rule_controller( ...
        [-0.3 13 15 13 -0.3],healthyV,normalRise,zeroVrate);
    assert(q(3) == 0 && q(2) > 0 && q(2) < 1 && q(4) == q(2));
    assert(~mode.insufficientRise);

    slowRise = normalRise;
    slowRise([1,5]) = 0.1;
    [q,mode] = q4_online_rule_controller( ...
        [-0.3 13 15 13 -0.3],healthyV,slowRise,zeroVrate);
    assert(all(q(2:4) == 1) && mode.insufficientRise);

    riskyV = healthyV;
    riskyV(3) = 0.40;
    [q,mode] = q4_online_rule_controller( ...
        [-0.3 13 15 13 -0.3],riskyV,normalRise,zeroVrate);
    assert(all(q(2:4) == 1) && mode.voltageRisk(3));

    iceV = healthyV;
    iceV(1) = 0.5;
    fallingV = zeroVrate;
    fallingV(1) = -0.03;
    [q,mode] = q4_online_rule_controller( ...
        [-0.3 13 15 13 -0.3],iceV,normalRise,fallingV);
    assert(q(2) == 1 && mode.iceRisk(1));
    assert(all(q >= 0 & q <= 1));
    fprintf('test_q4_online_rule_controller PASS\n');
end
