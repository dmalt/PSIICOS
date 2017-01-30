% Plotting part;
opacity = 0.2;
linewidth = 2;
m_radius = 0.003;

% ----------------- Psciicos --------------- %
subplot(3,3,1)
light_handle = con_ps.Plot(opacity, linewidth, m_radius);
view(0,0);
camlight(light_handle, 'left');

subplot(3,3,2)
light_handle = con_ps.Plot(opacity, linewidth, m_radius);
camlight(light_handle, 'left');
title('DICS PSIICOS');

subplot(3,3,3)
light_handle = con_ps.Plot(opacity, linewidth, m_radius);
view(180,0)
camlight(light_handle, 'right');

% --------------------- GCS ----------------- %
subplot(3,3,4)
light_handle = con_dics_gcs.Plot(opacity, linewidth, m_radius);
view(0,0);
camlight(light_handle, 'left');


subplot(3,3,5)
con_dics_gcs.Plot(opacity, linewidth, m_radius);
title('DICS GCS');

subplot(3,3,6)
light_handle = con_dics_gcs.Plot(opacity, linewidth, m_radius);
view(180,0)
camlight(light_handle, 'right');

% ----------------- DICS --------------- %
subplot(3,3,7)
light_handle = con_dics.Plot(opacity, linewidth, m_radius);
view(0,0)
camlight(light_handle, 'left');

subplot(3,3,8)
light_handle = con_dics.Plot(opacity, linewidth, m_radius);
title('DICS');

subplot(3,3,9)
light_handle = con_dics.Plot(opacity, linewidth, m_radius);
view(180,0)
camlight(light_handle, 'right');
