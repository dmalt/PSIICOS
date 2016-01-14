function CT_restored = RestoreSensorDimension(CT, UP)
	CT_reshaped = reshape(CT, sqrt(size(CT, 1)), sqrt(size(CT,1)));
    CT_restored = UP' * CT_reshaped * UP;