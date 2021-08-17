#version 300 es

precision highp float;

layout(location = 0) out vec4 out_acc;
layout(location = 1) out vec4 out_mask;
layout(location = 2) out vec4 out_ro;
layout(location = 3) out vec4 out_rd;
layout(location = 4) out vec4 out_ratts;

uniform sampler2D u_acc;
uniform sampler2D u_mask;
uniform sampler2D u_ro;
uniform sampler2D u_rd;
uniform sampler2D u_ratts;
uniform sampler2D u_env;

uniform uint u_frame;
uniform float u_time;
uniform vec2 u_resolution;

#define PI          3.1415926535897932384626433832795f
#define TWO_PI			 6.283185307179586476925286766559f

//internal RNG state 
uvec4 s0, s1; 
ivec2 pixel;

void rng_initialize(vec2 p, uint frame)
{
    pixel = ivec2(p);

    //white noise seed
    s0 = uvec4(p, uint(frame), uint(p.x) + uint(p.y));
    
    //blue noise seed
    s1 = uvec4(frame, frame*15843u, frame*31u + 4566u, frame*2345u + 58585u);
}

// https://www.pcg-random.org/
uvec4 pcg4d(inout uvec4 v)
{
	  v = v * 1664525u + 1013904223u;
    v.x += v.y*v.w; v.y += v.z*v.x; v.z += v.x*v.y; v.w += v.y*v.z;
    v = v ^ (v>>16u);
    v.x += v.y*v.w; v.y += v.z*v.x; v.z += v.x*v.y; v.w += v.y*v.z;
    return v;
}

float rand(){ return float(pcg4d(s0).x)/float(0xffffffffu); }
vec2 rand2(){ return vec2(pcg4d(s0).xy)/float(0xffffffffu); }
vec3 rand3(){ return vec3(pcg4d(s0).xyz)/float(0xffffffffu); }
vec4 rand4(){ return vec4(pcg4d(s0))/float(0xffffffffu); }

const struct Ray { vec3 o, d; };

float disp( in vec3 p )
{
	return pow( 0.5 + 0.5*cos( 21.0*p.x  + 1.5)*
					            sin( 23.0*p.y  + 2.0 )*
					            sin( 19.0*p.z  + 1.0 ), 3.0);
}

float obj( in vec3 p )
{
	vec3 ax = vec3(-2.0,2.0,1.0)/3.0;
	vec3 ce = vec3(0.0,-0.2,-0.2);

	float d1 = dot(p,ax) - 0.1;
    float d2 = length(p) - 1.0;
	float d3 = length( p-ce - ax*dot(p-ce,ax)) - 1.0;

	return max( max( d1, d2 ), -d3 );
}

vec2 map( vec3 p ){
  float d1 = obj( p );
  float d2 = obj( p*vec3(-1.0,-1.0,1.0) );

  vec2        res = vec2( d1, 0.0 );
  if( d2<d1 ) res = vec2( d2, 1.0 );

  res.x -= 0.1*disp( p );

  return res;
}

vec3 calcNormal( in vec3 pos )
{
    const vec3 eps = vec3(5e-4,0.0,0.0);

    return normalize( vec3(
      map( pos+eps.xyy ).x - map( pos-eps.xyy ).x,
      map( pos+eps.yxy ).x - map( pos-eps.yxy ).x,
      map( pos+eps.yyx ).x - map( pos-eps.yyx ).x ) );
}

bool march(inout vec3 o, vec3 dir){
    vec3 p = o;
    float e = 0.0;
    while(true){
        float d = 0.4*map(p).x;
        p += d*dir;
        e += d;
        if(d < 1e-4 || e > 12.0)
            break;
    }
    o = p;
    return e<=12.0;
}

float G(float dotNV, float k){
	return 1.0/(dotNV*(1.0f-k)+k);
}

// from http://filmicworlds.com/blog/optimizing-ggx-shaders-with-dotlh/
float ggx(vec3 N, vec3 V, vec3 L, float roughness, float F0){
	float alpha = roughness*roughness;

	vec3 H = normalize(V+L);

	float dotNL = clamp(dot(N,L),0.,1.);
	float dotNV = clamp(dot(N,V),0.,1.);
	float dotNH = clamp(dot(N,H),0.,1.);
	float dotLH = clamp(dot(L,H),0.,1.);

	float F, D, vis;

	float alphaSqr = alpha*alpha;
	float pi = 3.14159;
	float denom = dotNH * dotNH *(alphaSqr - 1.0) + 1.0;
	D = alphaSqr/(pi * denom * denom);

	float dotLH5 = pow(1.0 - dotLH, 5.0);
	F = F0 + (1.0 - F0)*(dotLH5);

	float k = alpha * 0.5;

	return dotNL * D * F * G(dotNL,k)*G(dotNV,k);
}
vec3 uniformSphere(){
	float phi = rand()*TWO_PI;
	float z = rand()*2.0f - 1.0f;
	float r = sqrt(max(1.0f - z * z, 0.0f));

	return vec3(cos(phi)*r, sin(phi)*r, z);
}
vec4 radiance(
    inout Ray r,
    inout vec4 mask,
    inout float depth)
{
  const vec3 l = vec3(0, 0, 5);

  if(map(r.o).x > 0.0){
    if(!march( r.o, r.d )){
      depth = 0.0;
      return vec4(1);
    }

    vec3 n = calcNormal(r.o);

    vec3 ld = normalize(l-r.o);

    mask.w = exp(-distance(l, r.o))*100.0;
    mask.y = 0.5+0.5*dot(n, ld);
    mask.z = ggx(n, -r.d, ld, 0.3, 0.7);

    r.d = refract(r.d, n, 1.0f/1.68);
  }

  ++depth;

  {
    // float hl = -log(1.0 - rand(h)) / m_sigmaT[int(rand()*3.0)];
    float hl = rand()*0.1;

    mask.x *= exp(-hl*2.0);

    r.o = r.o + r.d*hl;
    vec2 h = map(r.o);
    if(h.x > 0.0){
      // float att = 1.0+pow(distance(l, ray.pos), 2.0);

      float ssha = max(0.0, dot(calcNormal(r.o), normalize(l-r.o)));
      float diff = mix(mask.y, mix(mask.x, ssha, 0.2), mask.x);
      vec3 albedo = mix(vec3(0.6, 0.2, 0.2), vec3(0.2, 0.2, 0.6), h.y);
      depth = 0.0;
      return vec4(mask.w*(albedo*diff+0.025*mask.z), 1.0);
    } else {
      r.d = uniformSphere();
    }
  }

  return vec4(0);
}

mat3 setCamera()
{
  vec3 cz = normalize(vec3(0,0,-1));
  vec3 cx = normalize(cross(vec3(0, 1, 0), cz));
  vec3 cy = cross(cz, cx);
  return mat3(cx, cy, cz);
}

void main()
{
    // vec2 fc = gl_FragCoord.xy/u_resolution.xy;
    ivec2 ifc = ivec2(gl_FragCoord.xy);
  
    vec4 o_mask = vec4(1.0);
    vec4 o_ro = texelFetch(u_ro, ifc, 0);

    Ray r = Ray(vec3(0, 0, 1.65), vec3(0,0,1));
    rng_initialize(gl_FragCoord.xy, u_frame);

    if (o_ro.w == 0.0) {
      vec2 of = rand2()*2.0 - 1.0;
      vec2 p = (-u_resolution.xy + 2.0 * (gl_FragCoord.xy + of)) / u_resolution.y;      
      r.d = setCamera() * normalize(vec3(vec2(.5 * 0.036) * p, 0.019));
    }
    else {  
        o_mask = texelFetch(u_mask, ifc, 0);

        r = Ray(o_ro.xyz, texelFetch(u_rd, ifc, 0).xyz);
    }

    out_acc = texelFetch(u_acc, ifc, 0) + radiance(r, o_mask, o_ro.w);
    out_mask = o_mask;
    out_ro = vec4(r.o, o_ro.w);
    out_rd.xyz = r.d;
}
