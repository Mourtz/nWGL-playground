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

uniform float u_time;
uniform vec2 u_resolution;

#define EPS 1e-4
#define PI          3.1415926535897932384626433832795f
#define TWO_PI			 6.283185307179586476925286766559f

float seed = 0.;
float rand() { return fract(sin(seed+=0.1) * 43758.5453123); }
float rand(float h) { return fract(sin(seed+=h) * 43758.5453123); }
vec2 rand2() { float r0 = rand(); return vec2(r0, fract(sin(seed+=r0) * 43758.5453123)); }

const struct Ray { vec3 o, d; };

vec3 uniformSphere(){
	float phi = rand()*TWO_PI;
	float z = rand()*2.0f - 1.0f;
	float r = sqrt(max(1.0f - z * z, 0.0f));
	return vec3(cos(phi)*r, sin(phi)*r, z);
}

// by @iq
float Mandelbulb( vec3 p ){
    vec3 w = p;
    float m = dot(w,w);

    vec4 trap = vec4(abs(w),m);
  float dz = 1.0;


  for( int i=0; i<4; ++i )
  {
    dz = 8.0*pow(m,3.5)*dz + 1.0;

    float r = length(w);
    float b = 8.0*acos( clamp(w.y/r, -1.0, 1.0));
    float a = 8.0*atan( w.x, w.z );
    w = p + pow(r,8.0) * vec3( sin(b)*sin(a), cos(b), sin(b)*cos(a) );

    trap = min( trap, vec4(abs(w),m) );

    m = dot(w,w);
    if( m > 4.0 )
      break;
  }
  trap.x = m;

  return 0.25*log(m)*sqrt(m)/dz;
}

float map( vec3 p ){ return Mandelbulb(p); }

vec3 calcNormal( in vec3 pos )
{
    const vec3 eps = vec3(EPS,0.0,0.0);

    return normalize( vec3(
      map( pos+eps.xyy ) - map( pos-eps.xyy ),
      map( pos+eps.yxy ) - map( pos-eps.yxy ),
      map( pos+eps.yyx ) - map( pos-eps.yyx ) ) );
}


float intersect( in vec3 ro, in vec3 rd )
{
    float res = -1.0;
    float tmax = 16.0;
    float t = 0.01;
    while(true)
    {
        float h = map(ro+rd*t);
        if( h<EPS || t>tmax ) break;
        t += h;
    }
    
    if( t<tmax ) res = t;

    return res;
}

const vec3 m_density = vec3(0.01, 0.95, 1.0);

vec3 radiance(
    inout Ray r,
    inout vec3 mask,
    inout float depth)
{
  ++depth;

  vec3 acc = vec3(0);

  if(map(r.o) > 0.0){
    float t = intersect( r.o, r.d );

    if(t < 0.0){
      acc += mask*float(depth > 1.0)*2.0;
      depth = 0.0;
      return acc;
    }

    r.o = r.o + r.d*t;
    r.d = refract(r.d, calcNormal(r.o), 1.0/1.3);
  }

  const vec3 m_sigmaA = 0.6*m_density;
  const vec3 m_sigmaS = 1.2*m_density;
  const vec3 m_sigmaT = m_sigmaA+m_sigmaS;

  {
    float h = map(r.o);

    // float hl = -log(1.0 - rand(h)) / m_sigmaT[int(rand()*3.0)];
    float hl = abs(h-(0.7-rand()))*0.5;

    r.o = r.o + r.d*hl;

    h = map(r.o);
    mask *= exp(-min(hl, abs(h))*m_sigmaT);

    if(h > 0.0){

    } else {
      r.d = refract(r.d, vec3(0,sign(r.d.y),1), 1.0/1.5);
      mask *= m_sigmaS;
    }
  }

  float roulettePdf = max(max(mask.x, mask.y), mask.z);
  if (roulettePdf < 0.1f && depth > 2.0) {
      if (rand(depth*12.762) < roulettePdf)
          mask /= roulettePdf;
      else
          depth = 0.0;
  }

  return acc;
}

mat3 setCamera()
{
  vec3 cz = normalize(vec3(0,0,1));
  vec3 cx = normalize(cross(vec3(0, 1, 0), cz));
  vec3 cy = cross(cz, cx);
  return mat3(cx, cy, cz);
}

void main()
{
    // vec2 fc = gl_FragCoord.xy/u_resolution.xy;
    ivec2 ifc = ivec2(gl_FragCoord.xy);
  
    vec3 o_mask = vec3(1.0);
    vec3 o_ratts = texelFetch(u_ratts, ifc, 0).xyz;

    //seed = u_time + u_resolution.y * gl_FragCoord.x / u_resolution.x + gl_FragCoord.y / u_resolution.y;
    seed = o_ratts.x == 0.0 ? dot( gl_FragCoord.xy, vec2(12.9898, 78.233) ) + 1113.1 : o_ratts.y;
    seed *= sin(66.666*u_time);

    Ray r = Ray(vec3(0, 0, -1.65), vec3(0,0,1));

    if (o_ratts.x == 0.0) {
      vec2 of = rand2()*2.0 - 1.0;
      vec2 p = (-u_resolution.xy + 2.0 * (gl_FragCoord.xy + of)) / u_resolution.y;      
      r.d = setCamera() * normalize(vec3(vec2(.5 * 0.036) * p, 0.019));
    }
    else {  
        o_mask = texelFetch(u_mask, ifc, 0).xyz;

        r = Ray(texelFetch(u_ro, ifc, 0).xyz, texelFetch(u_rd, ifc, 0).xyz);
    }

    out_acc.xyz = texelFetch(u_acc, ifc, 0).xyz + radiance(r, o_mask, o_ratts.x);
    out_mask.xyz = o_mask;
    out_ro.xyz = r.o;
    out_rd.xyz = r.d;
    out_ratts.xy = vec2(o_ratts.x, seed);
}
