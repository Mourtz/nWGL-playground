#version 300 es

precision highp float;

layout(location = 0) out vec4 out_acc;
layout(location = 1) out vec4 out_mask;
layout(location = 2) out vec4 out_ro;
layout(location = 3) out vec4 out_rd;

uniform sampler2D u_acc;
uniform sampler2D u_mask;
uniform sampler2D u_ro;
uniform sampler2D u_rd;
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
        float h = 0.8*map(ro+rd*t);
        if( h<EPS || t>tmax ) break;
        t += h;
    }
    
    if( t<tmax ) res = t;

    return res;
}

float schlick(Ray r, vec3 n, float nc, float nt){
  float R0 = pow((nc - nt) / (nc + nt), 2.);
  return R0 + (1. - R0) * pow(1. + dot(n, r.d), 5.);
}

vec4 radiance(
    inout Ray r,
    inout vec3 mask,
    inout float depth)
{
  ++depth;

  const float nt = 2.418;
  vec4 res = vec4(0,0,0,0);

  float dist = map(r.o);
  if(dist > 0.0){
    float t = intersect( r.o, r.d );

    if(t < 0.0){
      depth = 0.0;
      return vec4(mask,1.0);
    }

    r.o = r.o + r.d*t;
    r.d = normalize(refract(r.d, calcNormal(r.o), 1.0f/nt));
    // r.d = uniformSphere();
    mask*=vec3(0.8, 0.7, 1.0);
    res = vec4(mask, 1.0);
  }
  else {
    // float hl = -log(1.0 - rand(h)) / m_sigmaT[int(rand()*3.0)];
    float hl = (1.0-rand())*0.2;
    mask *= exp(-hl*vec3(1.0,0.2,0.1)*20.0);

    float b = dot(mask,vec3(0.33));

    r.o = r.o + r.d*hl;
    if(abs(hl) >= -dist){
      res += vec4(mask, 0.2+0.8*(1.0-b));
    }else{
      r.d = uniformSphere();
      mask*=1.1;
      // r.d = normalize(refract(r.d, vec3(0,0,-1), 1.0f/nt));
    }

    //russian roulette
    if (b < 0.05f) {
      if (rand() < b){
        mask /= b;
      } else {
        depth = 0.0;
      }
    }

  }

  return res;
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

    vec4 ro = texelFetch(u_ro, ifc, 0);
    vec4 rd = texelFetch(u_rd, ifc, 0);

    //seed = u_time + u_resolution.y * gl_FragCoord.x / u_resolution.x + gl_FragCoord.y / u_resolution.y;
    seed = ro.w == 0.0 ? dot( gl_FragCoord.xy, vec2(12.9898, 78.233) ) + 1113.1 : rd.w;
    seed *= sin(66.666*u_time);

    Ray r = Ray(vec3(0, 0, -1.65), vec3(0,0,1));

    if (ro.w == 0.0) {
      vec2 of = rand2()*2.0 - 1.0;
      vec2 p = (-u_resolution.xy + 2.0 * (gl_FragCoord.xy + of)) / u_resolution.y;      
      r.d = setCamera() * normalize(vec3(vec2(.5 * 0.036) * p, 0.019));
    }
    else {  
        o_mask = texelFetch(u_mask, ifc, 0).xyz;

        r = Ray(ro.xyz, rd.xyz);
    }

    out_acc = texelFetch(u_acc, ifc, 0) + radiance(r, o_mask, ro.w);
    out_mask.xyz = o_mask;
    out_ro = vec4(r.o, ro.w);
    out_rd = vec4(r.d, seed);
}
