function Wmat = map_on_default(srcSurfMat, destSurfMat)
% ------------------------------------------------------------
% Generate linear operator to map data from individual cortex
% to average brain
% ------------------------------------------------------------
    nSrc  = size(srcSurfMat.Vertices, 1);
    nDest = size(destSurfMat.Vertices, 1);

    % return;
    nbNeighbors = 8;
    % Allocate interpolation matrix
    Wmat = spalloc(nDest, nSrc, nbNeighbors * nDest);
    % Split hemispheres
    [rHsrc, lHsrc, isConnected(1)]  = tess_hemisplit(srcSurfMat);
    [rHdest, lHdest, isConnected(2)] = tess_hemisplit(destSurfMat);
    % Get vertices
    srcVert  = double(srcSurfMat.Reg.Sphere.Vertices);
    destVert = double(destSurfMat.Reg.Sphere.Vertices);
    % If hemispheres are connected: process all at once
    if any(isConnected)
        rHsrc  = 1:nSrc;
        rHdest = 1:nDest;
        lHsrc  = [];
        lHdest = [];
    end
    % Re-interpolate using the sphere and the shepards algorithm
    WmatTmp = bst_shepards(destVert(rHdest,:), srcVert(rHsrc,:), nbNeighbors, 0);
    Wmat(rHdest,rHsrc) = WmatTmp;
    WmatTmp = bst_shepards(destVert(lHdest,:), srcVert(lHsrc,:), nbNeighbors, 0);
    Wmat(lHdest,lHsrc) = WmatTmp;
end
