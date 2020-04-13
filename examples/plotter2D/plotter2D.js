let stats = new Stats();
stats.showPanel(0);
document.body.appendChild(stats.dom);

let sandbox = new nWGL.main({ "fullscreen": true });

//------------------------- Textures -------------------------
// sandbox.addTexture({ "url": "tears_of_steel_bridge_8k.png" }, "hdr");

sandbox.addTexture({ "internalformat": "RGBA32F" }, "acc");
sandbox.addTexture({ "internalformat": "RGBA32F" }, "mask");

sandbox.addTexture({ "internalformat": "RGBA32F" }, "ro");
sandbox.addTexture({ "internalformat": "RGBA32F" }, "rd");


//------------------------- Shaders -------------------------
sandbox.addShader("../../shaders/vert.glsl", "vertex_shader", true);
sandbox.addShader("plotter2D.glsl", "plotter");

//------------------------- Plotting Program -------------------------
sandbox.addProgram(["vertex_shader", "plotter"], "plotter");

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

