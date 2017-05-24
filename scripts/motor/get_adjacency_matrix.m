function adj_mat = get_adjacency_matrix(Ctx)
% ------------------------------------------
% Get adjacency matrix for cortex surface
% ------------------------------------------
    verts = Ctx.Vertices;
    faces = Ctx.Faces;

    n_verts = length(verts);
    n_faces = length(faces);

    adj_mat = sparse(n_verts, n_verts);
    for i_face = 1:n_faces
        face = faces(i_face,:);
        if ~adj_mat(face(1),face(2))
            adj_mat(face(1),face(2)) = 1;
            adj_mat(face(2),face(1)) = 1;
        end

        if ~adj_mat(face(1), face(3))
            adj_mat(face(1), face(3)) = 1;
            adj_mat(face(3), face(1)) = 1;
        end

        if ~adj_mat(face(2), face(3))
            adj_mat(face(2), face(3)) = 1;
            adj_mat(face(3), face(2)) = 1;
        end
    end
end
