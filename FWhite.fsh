#ifdef GL_ES
// Set default precision to medium
precision mediump int;
precision mediump float;
#endif

varying vec3 vc_color;

void main()
{
    gl_FragColor = vec4(vc_color, 1.0);
}
