#version 300 es

precision highp float;

layout(location = 0) in vec2 a_position;

void main(void) {
    gl_Position.xy = a_position;
}
