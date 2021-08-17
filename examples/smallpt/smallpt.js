let stats = new Stats();
stats.showPanel(0);
document.body.appendChild(stats.dom);

let sandbox = new nWGL.main({ "fullscreen": true });

//------------------------- Framebuffers -------------------------
let fb = sandbox.addFrameBuffer({ "internalformat": "RGBA32F", "totalBuffers": 4 }, "backbuffer");

//------------------------- Textures -------------------------
// sandbox.addTexture({ "url": "tears_of_steel_bridge_8k.png" }, "hdr");

sandbox.addTexture({ "internalformat": "RGBA32F" }, "acc");
sandbox.addTexture({ "internalformat": "RGBA32F" }, "mask");

sandbox.addTexture({ "internalformat": "RGBA32F" }, "ro");
sandbox.addTexture({ "internalformat": "RGBA32F" }, "rd");


//------------------------- Shaders -------------------------
sandbox.addShader("../../shaders/vert.glsl", "vertex_shader", true);
sandbox.addShader("smallpt.glsl", "pt_shader");
sandbox.addShader("../../shaders/display.glsl", "display_shader");

//------------------------- RayTracing Program -------------------------
sandbox.addProgram(["vertex_shader", "pt_shader"], "raytracing");

sandbox.programs["raytracing"].addUniform("u_acc", "1i", 0);
sandbox.programs["raytracing"].addUniform("u_mask", "1i", 1);
sandbox.programs["raytracing"].addUniform("u_ro", "1i", 2);
sandbox.programs["raytracing"].addUniform("u_rd", "1i", 3);
// sandbox.programs["raytracing"].addUniform("u_env", "1i", 5);

//------------------------- Display Program -------------------------
sandbox.addProgram(["vertex_shader", "display_shader"], "display");

// sandbox.programs["display"].addUniform("u_cont", "1f", 1);
sandbox.programs["display"].addUniform("u_tex", "1i", 0);

//------------------------- Render Passes -------------------------
sandbox.composer.addPass(
  {
    "render": function () {
      sandbox.program = "raytracing";

      sandbox.setTexture("u_acc", sandbox.textures["acc"].tex, 0);
      sandbox.setTexture("u_mask", sandbox.textures["mask"].tex, 1);
      sandbox.setTexture("u_ro", sandbox.textures["ro"].tex, 2);
      sandbox.setTexture("u_rd", sandbox.textures["rd"].tex, 3);

      sandbox.bindFramebuffer(fb.fb);
      fb.setTexture(0, fb.t0);
      fb.setTexture(1, fb.t1);
      fb.setTexture(2, fb.t2);
      fb.setTexture(3, fb.t3);
    },
    "swapBuffer": false
  },
  {
    "render": function () {
      sandbox.program = "display";
      // sandbox.uniform("u_cont", 1.0 / (sandbox.frame + 1));
      sandbox.setTexture("u_tex", fb.t0, 0);
      sandbox.bindFramebuffer(null);
    }
  },
  {
    "compute": function () {
      // texture ping pong
      sandbox.textures["acc"].swap(fb.textures[0]);
      sandbox.textures["mask"].swap(fb.textures[1]);
      sandbox.textures["ro"].swap(fb.textures[2]);
      sandbox.textures["rd"].swap(fb.textures[3]);
    }
  }
);

//------------------------- Render Loop -------------------------
(function render() {
  stats.begin();
  sandbox.draw();
  stats.end();
  window.requestAnimationFrame(render);
})();

// if(false){
//   render();
// }else{
//   // halt until the texture loads
//   {
//     let t = performance.now();
//     (function check(){
//       if(performance.now() - t > 500){
//         if(sandbox.textures["hdr"].image.complete){
//           sandbox.program = "raytracing";
//           sandbox.setTexture("u_env", sandbox.textures["hdr"].tex, 5);
    
//           render();
//           return;
//         }
//         t = performance.now();
//       }
//       window.requestAnimationFrame(check);
//     })();
//   }
// }

