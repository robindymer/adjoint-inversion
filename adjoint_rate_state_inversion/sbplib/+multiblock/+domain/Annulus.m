classdef Annulus < multiblock.DefCurvilinear
    properties
        r_inner % Radii of inner disk
        c_inner % Center of inner disk
        r_outer % Radii of outer disk
        c_outer % Radii of outer disk
    end

    methods
        function obj = Annulus(r_inner, c_inner, r_outer, c_outer)
            default_arg('r_inner', 0.3);
            default_arg('c_inner', [0; 0]);
            default_arg('r_outer', 1)
            default_arg('c_outer', [0; 0]);
            % Assert that the problem is well-defined
            d = norm(c_outer-c_inner,2);
            assert(r_inner > 0, 'Inner radius must be greater than zero');
            assert(r_outer > d+r_inner, 'Inner disk not contained in outer disk');
            
            cir_out_A = parametrization.Curve.circle(c_outer,r_outer,[-pi/2 pi/2]);
            cir_in_A = parametrization.Curve.circle(c_inner,r_inner,[pi/2 -pi/2]);
            
            cir_out_B = parametrization.Curve.circle(c_outer,r_outer,[pi/2 3*pi/2]);
            cir_in_B = parametrization.Curve.circle(c_inner,r_inner,[3*pi/2 pi/2]);

            c0_out = cir_out_A(0);
            c1_out = cir_out_A(1);
            
            c0_in_A = cir_in_A(1);
            c1_in_A = cir_in_A(0);
            
            c0_out_B = cir_out_B(0);
            c1_out_B = cir_out_B(1);
            
            c0_in_B = cir_in_B(1);
            c1_in_B = cir_in_B(0);


            sp2_A = parametrization.Curve.line(c0_in_A,c0_out);
            sp3_A = parametrization.Curve.line(c1_in_A,c1_out);
            
            sp2_B = parametrization.Curve.line(c0_in_B,c0_out_B);
            sp3_B = parametrization.Curve.line(c1_in_B,c1_out_B);


            A = parametrization.Ti(sp2_A, cir_out_A, sp3_A.reverse, cir_in_A); 
            B = parametrization.Ti(sp2_B , cir_out_B,sp3_B.reverse, cir_in_B );
            
            blocks = {A,B};
            blocksNames = {'A','B'};

            conn = cell(2,2);

            conn{1,2} = {'n','s'};
            conn{2,1} = {'n','s'};

            boundaryGroups = struct();
            boundaryGroups.out = multiblock.BoundaryGroup({{1,'e'},{2,'e'}});
            boundaryGroups.in = multiblock.BoundaryGroup({{1,'w'},{2,'w'}});
            boundaryGroups.all = multiblock.BoundaryGroup({{1,'e'},{2,'w'},{1,'w'},{2,'e'}});

            obj = obj@multiblock.DefCurvilinear(blocks, conn, boundaryGroups, blocksNames);

            obj.r_inner = r_inner;
            obj.r_outer = r_outer;
            obj.c_inner = c_inner;
            obj.c_outer = c_outer;
        end

        function ms = getGridSizes(obj, m)
            mx = m;
            % Use same grid spacing along inner 
            % half circle as in radial direction
            ds = pi*(obj.r_inner);
            dr = (obj.r_outer-obj.r_inner);
            my = ceil(ds/dr*(mx-1))+1;

            ms = {[mx my], [mx my]};
        end
    end
end
