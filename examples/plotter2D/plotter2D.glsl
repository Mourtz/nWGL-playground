#version 300 es

precision highp float;

layout(location = 0) out vec4 out_color;

uniform float u_time;
uniform vec2 u_resolution;
uniform vec3 u_mouse;

#define time u_time*0.001

#define AA 2	// make this 1 if you have a slow computer
#define MOUSE_SENSITIVITY vec2(1, 1)

#define PI 3.141592653589793238

#define dStep 1e-3
#define deriv(func, x, y0) (func(x+dStep)-y0)/dStep

vec3 calc(in vec2 fragCoord, in vec2 mouse, in float px)
{  
	float x, y, dy, ddy;
    
#define ZOOM 7.0 // 10-30
    
    vec2 uv = (mouse + (fragCoord/u_resolution.xy-0.5)) * ZOOM;    
    
    // grid & axis
    vec2 cb = step(vec2(ZOOM*px), fract(uv)) * step(ZOOM*(px*2.0), abs(uv));
    vec3 col = vec3(cb.x*cb.y);
            
    // 1st function(green)
	  x = uv.x + time;
    y = cos(x);
    dy = deriv(cos, x, y);
    y -= uv.y;
    col = mix(vec3(0.0, 1.0, 0.0), col, smoothstep(0.0, ZOOM*(px*2.0), abs(y) / (abs(dy) + 1.0)));
    // 1st derivative (yellow)
    dy -= uv.y;
    col = mix(vec3(1.0, 1.0, 0.0), col, smoothstep(0.0, ZOOM*(px*2.0), abs(dy) / (dFdx(dy) + 1.8)));

    // 2nd function(red)
    x = -uv.x*sin(time);
    y = exp(x);
    dy = deriv(exp, x, y);
    y = 3.0 + y - uv.y;
    col = mix(vec3(1.0, 0.0, 0.0), col, smoothstep(0.0, ZOOM*(px*2.0), abs(y) / (abs(dy) + 1.0) ));
    
    // 3rd function(blue)
    x = time+sin(uv.x);
    y = sin(x);
    dy = deriv(sin, x, y);
    y = -4.0 + y - uv.y;
    col = mix(vec3(0.0, 0.0, 1.0), col, smoothstep(0.0, ZOOM*(px*2.0), abs(y) / (abs(dy) + 1.0)));
    // 1st derivative (orange)
    dy = -4.0 + dy - uv.y;
    col = mix(vec3(1.0, 0.65, 0.0), col, smoothstep(0.0, ZOOM*(px*2.0), abs(dy) / (dFdx(dy) + 1.8)));
    
#undef ZOOM
    
    return col;
}

void main()
{  
    vec2 fragCoord = gl_FragCoord.xy;
    fragCoord.x*=u_resolution.x/u_resolution.y;
    vec2 mouse = (u_mouse.xy-0.5)*MOUSE_SENSITIVITY;
    
    vec3 col = vec3(0.0);
#if AA > 1
    const float aa = float(AA);
    
    for( float j=0.0; j<aa; ++j )
	for( float i=0.0; i<aa; ++i )
	{
		vec2 of = -0.5 + vec2(i, j)/aa;
	    col += calc(fragCoord+of, mouse, 1.0/u_resolution.y);
	}
	col /= (aa*aa);
#else
	col = calc(fragCoord, mouse, 1.0/u_resolution.y);
#endif
    
    out_color = vec4(col, 1.0);
}