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

uniform float u_time;
uniform vec2 u_resolution;

#define EPS 1e-3
#define PI 3.14159265359
#define DIFF 0.0
#define SPEC 1.0
#define REFR 2.0
#define NUM_SPHERES 9

float seed = 0.;
float rand() { return fract(sin(seed+=0.1) * 43758.5453123); }
vec2 rand2() { float r0 = rand(); return vec2(r0, fract(sin(seed+=r0) * 43758.5453123)); }

const struct Ray {
    vec3 o, d;
};
const struct Sphere {
    float r;
    vec3 p, e, c;
    float refl;
};

const Sphere lightSourceVolume = Sphere(20., vec3(50., 81.6, 81.6), vec3(12.), vec3(0.), DIFF);

#define S_R 1e3
const Sphere spheres[NUM_SPHERES] = Sphere[](
    Sphere(S_R, vec3(-S_R + 1., 40.8, 81.6), vec3(0.), vec3(.75, .25, .25), DIFF),
    Sphere(S_R, vec3(S_R + 99., 40.8, 81.6), vec3(0.), vec3(.25, .25, .75), DIFF),
    Sphere(S_R, vec3(50., 40.8, -S_R), vec3(0.), vec3(.75), DIFF),
    Sphere(S_R, vec3(50., 40.8, S_R + 170.), vec3(0.), vec3(0.75), DIFF),
    Sphere(S_R, vec3(50., -S_R, 81.6), vec3(0.), vec3(.75), DIFF),
    Sphere(S_R, vec3(50., S_R + 81.6, 81.6), vec3(0.), vec3(.75), DIFF),
    Sphere(16.5, vec3(27., 16.5, 47.), vec3(0.), vec3(1.), SPEC),
    Sphere(16.5, vec3(73., 16.5, 78.), vec3(0.), vec3(1.), REFR),
    Sphere(600., vec3(50., 681.33, 81.6), vec3(12.), vec3(0.), DIFF));

float intersect(Sphere s, Ray r)
{
    vec3 op = s.p - r.o;
    float t, b = dot(op, r.d), det = b * b - dot(op, op) + s.r * s.r;
    if (det < 0.)
        return 0.;
    else
        det = sqrt(det);
    return (t = b - det) > EPS ? t : ((t = b + det) > EPS ? t : 0.);
}

int intersect(Ray r, out float t, out Sphere s)
{
    int id = -1;
    t = 1e5;
    s = spheres[0];
    for (int i = 0; i < NUM_SPHERES; ++i) {
        Sphere S = spheres[i];
        float d = intersect(S, r);
        if (d != 0. && d < t) {
            t = d;
            id = i;
            s = S;
        }
    }
    return id;
}

vec3 jitter(vec3 d, float phi, float sina, float cosa)
{
    vec3 w = normalize(d), u = normalize(cross(w.yzx, w)), v = cross(w, u);
    return (u * cos(phi) + v * sin(phi)) * sina + w * cosa;
}

vec3 radiance(
    inout Ray r,
    inout vec3 mask,
    inout float depth)
{
    vec3 acc = vec3(0.0);

    float t;
    Sphere obj;
    int id = intersect(r, t, obj);

    if (++depth > 64.0 || id < 0) {
        depth = 0.0;
        return acc;
    }

    vec3 x = t * r.d + r.o;
    vec3 n = normalize(x - obj.p), nl = n * sign(-dot(n, r.d));

    if (obj.refl == DIFF) {
        float r2 = rand();
        vec3 d = jitter(nl, 2. * PI * rand(), sqrt(r2), sqrt(1. - r2));
        vec3 e = vec3(0.);

        {
            Sphere s = lightSourceVolume;

            vec3 l0 = s.p - x;
            float cos_a_max = sqrt(1. - clamp(s.r * s.r / dot(l0, l0), 0., 1.));
            float cosa = mix(cos_a_max, 1., rand());
            vec3 l = jitter(l0, 2.*PI*rand(), sqrt(1. - cosa*cosa), cosa);

            if (intersect(Ray(x+nl*EPS, l), t, s) == 8) {
                float omega = 2. * PI * (1. - cos_a_max);
                e += (s.e * clamp(dot(l, n),0.,1.) * omega) / PI;
            }
        }
        acc += mask * obj.e + mask * obj.c * e;
        mask *= obj.c;
        r = Ray(x + nl * EPS, d);
    }
    else if (obj.refl == SPEC) {
        acc += mask * obj.e;
        mask *= obj.c;
        r = Ray(x + nl * EPS, reflect(r.d, n));
    }
    else {
        float a = dot(n, r.d), ddn = abs(a);
        float nc = 1., nt = 1.5, nnt = mix(nc / nt, nt / nc, float(a > 0.));
        float cos2t = 1. - nnt * nnt * (1. - ddn * ddn);
        r = Ray(x + nl * EPS, reflect(r.d, n));
        if (cos2t > 0.) {
            vec3 tdir = normalize(r.d * nnt + sign(a) * n * (ddn * nnt + sqrt(cos2t)));
            float R0 = (nt - nc) * (nt - nc) / ((nt + nc) * (nt + nc)),
                  c = 1. - mix(ddn, dot(tdir, n), float(a > 0.));
            float Re = R0 + (1. - R0) * c * c * c * c * c, P = .25 + .5 * Re, RP = Re / P, TP = (1. - Re) / (1. - P);
            if (rand() < P) {
                mask *= RP;
            }
            else {
                mask *= obj.c * TP;
                r = Ray(x - nl * EPS, tdir);
            }
        }
    }
  
    float roulettePdf = max(max(mask.x, mask.y), mask.z);
    if (roulettePdf < 0.1f && depth > 2.0) {
        if (rand() < roulettePdf)
            mask /= roulettePdf;
        else
            depth = 0.0;
    }
  
    return acc;
}

mat3 setCamera(in vec3 ro, in vec3 rt, in float cr)
{
    vec3 cw = normalize(rt - ro);
    vec3 cp = vec3(sin(cr), cos(cr), 0.0);
    vec3 cu = normalize(cross(cw, cp));
    vec3 cv = normalize(cross(cu, cw));
    return mat3(cu, cv, -cw);
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

    Ray r;

    if (o_ratts.x == 0.0) {
        vec2 of = -0.5 + rand2();
        vec2 p = (-u_resolution.xy + 2.0 * (gl_FragCoord.xy + of)) / u_resolution.y;
        mat3 ca = setCamera(r.o, vec3(0.0, 0.0, -1.5), 0.0);
        r = Ray(vec3(50.0, 40.8, 169.), normalize(ca * vec3(p, -1.5)));
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
