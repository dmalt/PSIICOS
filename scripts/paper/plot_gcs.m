% Plotting part;
opacity = 0.2;
linewidth = 2;
m_radius = 0.003;

subplot(2,3,1)
light_handle = con_ps.Plot(opacity, linewidth, m_radius);
view(0,0);
camlight(light_handle, 'left');

subplot(2,3,2)
light_handle = con_ps.Plot(opacity, linewidth, m_radius);
camlight(light_handle, 'left');
title('PSIICOS');

subplot(2,3,3)
light_handle = con_ps.Plot(opacity, linewidth, m_radius);
view(180,0)
camlight(light_handle, 'right');

subplot(2,3,4)
light_handle = con_dics_gcs.Plot(opacity, linewidth, m_radius);
view(0,0);
camlight(light_handle, 'left');


subplot(2,3,5)
con_dics_gcs.Plot(opacity, linewidth, m_radius);
title('DICS GCS');

subplot(2,3,6)
light_handle = con_dics_gcs.Plot(opacity, linewidth, m_radius);
view(180,0)
camlight(light_handle, 'right');




