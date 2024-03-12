function J = glut1(gGLUT,p)
% This function is going to be set up such that it is measured and never to
% be touched again. The measurement will come from the Lab. 

J = gGLUT * ((p.Tec_glu * p.glu_e ) / (p.Mec_glu + p.glu_e ));
end