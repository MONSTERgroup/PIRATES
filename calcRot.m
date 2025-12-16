function [rot_mat] = calcRot(it,pos)

p1i = vector3d([pos(it-1, 1, 1),  pos(it-1, 1, 2),  pos(it-1, 1, 3)]);
        p1f = vector3d([pos(it,   1, 1),  pos(it,   1, 2),  pos(it,   1, 3)]);
        p2i = vector3d([pos(it-1, 2, 1),  pos(it-1, 2, 2),  pos(it-1, 2, 3)]);
        p2f = vector3d([pos(it,   2, 1),  pos(it,   2, 2),  pos(it,   2, 3)]);
        p3i = vector3d([pos(it-1, 3, 1),  pos(it-1, 3, 2),  pos(it-1, 3, 3)]);
        p3f = vector3d([pos(it,   3, 1),  pos(it,   3, 2),  pos(it,   3, 3)]);

        %Set Point 1 to be the origin
        u1 = p2i - p1i;
        v1 = p2f - p1f;
        u2 = p3i - p1i;
        v2 = p3f - p1f;

        u1n = normalize(u1);
        v1n = normalize(v1);
        u2n = normalize(u2);
        v2n = normalize(v2);

        delta = abs(angle(u1n,u2n,'noSymmetry') - angle(v1n,v2n,'noSymmetry'));
        difference = max(delta(:))./degree;

        if delta ~= 0
            warning(['Inconsitent pairs of vectors! sThe angle between u1, u2 and v1, v2 needs ' ....
                'to be the same, but differs by ' num2str(difference),mtexdegchar, '\n']);
        end
        
        %Calculate Rotation Matrix from rotation of two axes
        %See mtex-toolbox.github.io/RotationDefinition.html for confirmation
        %that this method is valid.

        
            rot_mat = rotation.map(u1,v1,u2,v2);
            rot_mat = rot_mat.matrix;

            for ii = 1:length(rot_mat)
                for jj = 1:length(rot_mat)
                    if abs(rot_mat(ii,jj)) < 1e-3
                        rot_mat(ii,jj) = 0;
                    end
                end
            end
        