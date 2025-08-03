class Particle {
  constructor(z, y, x = 0) {
    this.position = { z, y, x };
  }
}

function configureSlider(params) {
  const slider = document.getElementById(params.sliderName);
  const sliderValue = document.getElementById(`${params.sliderName}Value`);
  params.setterMethod(params.defaultVal);
  sliderValue.textContent = Number(params.defaultVal).toFixed(2);
  slider.value = (100 / (params.maxVal - params.minVal)) * params.defaultVal;
  slider.addEventListener("input", () => {
    const val =
      params.minVal + ((params.maxVal - params.minVal) / 100) * slider.value;
    sliderValue.textContent = Number(val).toFixed(2);
    params.setterMethod(val);
  });
}

function makeGradient(x, xMin, xMax) {
  const t = Math.max(0, Math.min(1, (x - xMin) / (xMax - xMin))); // clamp to [0,1]

  let r, g, b;

  if (t < 0.5) {
    // pastel blue → pastel white
    const t2 = t * 2;
    r = Math.round(200 + 55 * t2); // 200–255
    g = Math.round(200 + 55 * t2); // 200–255
    b = 255;
  } else {
    // pastel white → pastel red
    const t2 = (t - 0.5) * 2;
    r = 255;
    g = Math.round(255 - 55 * t2); // 255–200
    b = Math.round(255 - 55 * t2); // 255–200
  }

  return [r, g, b];
}

class Rendered2D {
  constructor(sim) {
    this.sim = sim;
    this.canvas = document.getElementById("simulationCanvas");
    this.ctx = this.canvas.getContext("2d");
    this.domain_limits = sim.get_domain_limits();
    this.target_dt = 1 / 60;

    // Sliders
    configureSlider({
      sliderName: "viscositySlider",
      minVal: 0.01,
      maxVal: 0.2,
      defaultVal: 0.1,
      setterMethod: (visc) => this.sim.set_viscosity(visc),
    });
    configureSlider({
      sliderName: "gravitySlider",
      minVal: 0,
      maxVal: 20,
      defaultVal: 9.8,
      setterMethod: (g) => this.sim.set_gravity(g),
    });
    configureSlider({
      sliderName: "specificVolumeSlider",
      minVal: 1 / 100,
      maxVal: 1 / 10,
      defaultVal: 1 / 50,
      setterMethod: (k) => this.sim.set_specific_volume(1 / k),
    });
    this.densityMin = 500;
    this.densityMax = 2000;
    this.newParticlesDensity = 1000;
    configureSlider({
      sliderName: "densitySlider",
      minVal: this.densityMin,
      maxVal: this.densityMax,
      defaultVal: this.newParticlesDensity,
      setterMethod: (rho) => {
        this.newParticlesDensity = rho;
      },
    });

    const mouseInteractionSelector = document.getElementById(
      "mouseInteractionSelector"
    );
    this.mouseInteractionSelectorState = mouseInteractionSelector.value;
    mouseInteractionSelector.addEventListener("change", () => {
      this.mouseInteractionSelectorState = mouseInteractionSelector.value;
    });

    this.numberOfParticlesElement = document.getElementById(
      "numberOfParticlesValue"
    );
    this.numberOfParticlesElement.textContent = this.sim.get_particles().size();

    // Mouse events
    // I want it to be screenHeight / 6. So it's screenHeight / 6 / screenHeight * simY.
    this.isDragging = false;
    this.simDragRadius = 0.2 * this.domain_limits.y;
    this.lastX = 0;
    this.lastY = 0;
    this.dragButton = 0;
    this.canvas.addEventListener("mousedown", (e) => {
      this.isDragging = true;
      [this.lastX, this.lastY] = [e.offsetX, e.offsetY];
      this.dragButton = e.button;
    });

    this.canvas.addEventListener("mousemove", (e) => {
      if (!this.isDragging) return;
      [this.lastX, this.lastY] = [e.offsetX, e.offsetY];
    });

    this.canvas.addEventListener("mouseup", () => {
      this.isDragging = false;
    });

    this.canvas.addEventListener("mouseleave", () => {
      this.isDragging = false;
    });
  }

  width = () => {
    return this.canvas.width;
  };
  height = () => {
    return this.canvas.height;
  };

  applyCentralForce = () => {
    const magnitude = this.dragButton === 0 ? 30 : -30;
    const simCoords = this.coordinatesScreenToSim({
      z: 0,
      y: this.lastY,
      x: this.lastX,
    });
    this.sim.apply_central_force(simCoords, magnitude, this.simDragRadius);
  };

  addOrRemoveParticle = () => {
    const simCoords = this.coordinatesScreenToSim({
      z: 0,
      y: this.lastY,
      x: this.lastX,
    });
    if (this.dragButton === 0) {
      this.sim.add_particle(simCoords, this.newParticlesDensity);
    } else {
      this.sim.remove_particle_at(simCoords);
    }
    this.numberOfParticlesElement.textContent = this.sim.get_particles().size();
  };

  coordinatesSimToScreen = (simCoords) => {
    return {
      z: 0,
      y:
        ((this.domain_limits.y - simCoords.y) * this.height()) /
        this.domain_limits.y,
      x: (simCoords.x * this.width()) / this.domain_limits.x,
    };
  };

  coordinatesScreenToSim = (screenCoords) => {
    return {
      z: this.domain_limits.z / 2,
      y: (1 - screenCoords.y / this.height()) * this.domain_limits.y,
      x: (screenCoords.x / this.width()) * this.domain_limits.x,
    };
  };

  update = () => {
    this.sim.update(this.target_dt);
    return this.sim.get_particles();
  };

  draw = (particles) => {
    const ctx = this.ctx;
    ctx.fillStyle = "black";
    ctx.fillRect(0, 0, this.width(), this.height());

    for (let i = 0; i < particles.size(); i++) {
      const p = particles.get(i);

      const [r, g, b] = makeGradient(
        p.density,
        this.densityMin,
        this.densityMax
      );
      // ctx.fillStyle = `rgb(${color[0]},${color[1]},${color[2]})`;
      // ctx.fillStyle = `rgb(255,0,0)`;

      const pos = this.coordinatesSimToScreen(p.position);
      ctx.beginPath();
      ctx.arc(pos.x, pos.y, 4, 0, 2 * Math.PI);
      const radius = 4;
      const gradient = ctx.createRadialGradient(
        pos.x,
        pos.y,
        0,
        pos.x,
        pos.y,
        radius
      );
      gradient.addColorStop(0, `rgba(${r},${g},${b}, 1.0)`);
      gradient.addColorStop(1, `rgba(${r},${g},${b}, 0.7)`);

      ctx.fillStyle = gradient;
      ctx.beginPath();
      ctx.arc(pos.x, pos.y, radius, 0, 2 * Math.PI);
      ctx.fill();
    }
  };

  handleInteractions = () => {
    if (!this.isDragging) {
      return;
    }

    switch (this.mouseInteractionSelectorState) {
      case "Force":
        this.applyCentralForce();
        break;
      case "Particles":
        this.addOrRemoveParticle();
        break;
    }
  };

  loop = () => {
    this.handleInteractions();
    const particles = this.update();
    this.draw(particles);
    requestAnimationFrame(this.loop);
  };
}

SimModule().then((mod) => {
  document.getElementById("status").textContent = "Module loaded.";

  const sim = mod.initialize_simulation();
  const renderer = new Rendered2D(sim);
  renderer.loop();
});
