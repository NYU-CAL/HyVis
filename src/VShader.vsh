#ifdef GL_ES
// Set default precision to medium
precision mediump int;
precision mediump float;
#endif

uniform mat4 ortho_matrix;
uniform vec3 const_color;

attribute vec4 a_position;
attribute vec3 a_color;

varying vec3 v_color;
varying vec3 vc_color;

void main()
{
    gl_Position = ortho_matrix * a_position;
    v_color = a_color;
    vc_color = const_color;
}
