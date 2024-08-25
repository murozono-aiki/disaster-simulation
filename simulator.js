/**
 * 加速度を計算するときに使う項
 * @type {{pressureTerm:{x:number, y:number, z:number}, viscosityTerm:{x:number, y:number, z:number}}[]}
 */
const terms = [];

const h = 0.012; //影響半径
const particleMass = 0.0002; //粒子の質量

const g = createVector3(0, -9.8, 0);  // 重力加速度
const pressureStiffness = 200; //圧力係数
const restDensity = 1000; //静止密度
const viscosity = 1;  // 粘性係数

const densityCoef = particleMass * 315 / (64 * Math.PI * Math.pow(h,9)); //密度計算で使うヤツ

const pressureCoef = particleMass * 45 / (Math.PI * Math.pow(h,6)); //圧力項計算で使うヤツ
const viscosityCoef = viscosity * particleMass * 45 / (Math.PI * Math.pow(h,6)); //粘性項計算で使うヤツ

//壁は立方体
const wallWidth = 10; 
const thickness = 0.5/* 0.1 */; //壁となる粒子の厚み（個数）
const interval = h / thickness;


const shaft = ["x","y","z"];

//壁となる、速度がゼロで固定の粒子群
function makeWall(particles){
    let newPosition = createVector3(0,0,0);
    for(let x = 0; x < 11 ; x = x + 10){
        newPosition[shaft[2]] = x;
        for(let i=0 ; i< wallWidth / interval; i++){//高さ
            for(let j = 0; j < wallWidth / interval; j++){//横幅
                newPosition[shaft[1]] = (interval / 2) + i * interval;
                newPosition[shaft[0]] = (interval / 2) + j * interval;
                const newParticle = createParticle(newPosition, true);
                particles.push(newParticle);
            }

        }  
    }
    for(let x = 0; x < 11 ; x = x + 10){
        newPosition[shaft[0]] = x;
        for(let i=0 ; i< wallWidth / interval; i++){//高さ
            for(let j = 0; j < wallWidth / interval; j++){//横幅
                newPosition[shaft[1]] = (interval / 2) + i * interval;
                newPosition[shaft[2]] = (interval / 2) + j * interval;
                const newParticle = createParticle(newPosition, true);
                particles.push(newParticle);
            }

        }  
    }
    newPosition[shaft[1]] = 0;
        for(let i=0 ; i< wallWidth / interval; i++){//x方向
            for(let j = 0; j < wallWidth / interval; j++){//z方向
                newPosition[shaft[0]] = (interval / 2) + i * interval;
                newPosition[shaft[2]] = (interval / 2) + j * interval;
                const newParticle = createParticle(newPosition, true);
                particles.push(newParticle);
            }

        }
    //console.log(particles)
}
function addParticles(particles) {
    for (let i = 2; i < 5000; i += 2) {
        particles.push(createParticle(createVector3(1, i, 1)));
    }
}



/**
 * 粒子の密度計算
 * @param {{position:{x:number, y:number, z:number}, velocity:{x:number, y:number, z:number}, force:{x:number, y:number, z:number}, density:number, pressure:number}[]} particles - 粒子のリスト
 */
function calcDensity(particles) {
    const h2 = h*h; //事前にhの二乗を計算しておく
    for (let i = 0; i < particles.length; i++) { //一つづつ粒子の密度を計算

        reportProgress(reportHeader + "密度を計算中..." + i + "/" + particles.length);

        let nowParticle = particles[i]; //今回計算する粒子
        let sum = 0; //足し合わせる変数
        for (let j = 0; j < particles.length; j++) { //他の粒子全てについて
            if(i == j){continue;} //自分自身だったらスキップ
            let nearParticle = particles[j];
            
            let diff = subVector3(nearParticle.position, nowParticle.position); //粒子距離
            const r2 = dotVector3(diff, diff); //粒子距離の２乗

            //粒子距離がhより小さい場合だけ計算する
            if ( r2 < h2 ) {
                const c = h2 - r2;
                sum += Math.pow(c,3); //(h2-r2)の３乗
            }
        }

        nowParticle.density = sum * densityCoef; //密度が求まった
    }
}

/**
 * 粒子の圧力計算
 * @param {{position:{x:number, y:number, z:number}, velocity:{x:number, y:number, z:number}, force:{x:number, y:number, z:number}, density:number, pressure:number}[]} particles - 粒子のリスト
 */
function calcPressure(particles) {
    for(let i = 0; i < particles.length; i++) { //一つづつ粒子の圧力を計算

        reportProgress(reportHeader + "圧力を計算中..." + i + "/" + particles.length);

        particles[i].pressure = pressureStiffness * ( particles[i].density - restDensity );
    }
}

/**
 * 粒子の圧力項計算
 * @param {{position:{x:number, y:number, z:number}, velocity:{x:number, y:number, z:number}, force:{x:number, y:number, z:number}, density:number, pressure:number}[]} particles - 粒子のリスト
 */
function calcPressureTerm(particles) {
    const h2 = h*h; //事前にhの二乗を計算しておく
    for (let i = 0; i < particles.length; i++) { //一つづつ粒子の密度を計算

        reportProgress(reportHeader + "圧力項を計算中..." + i + "/" + particles.length);

        let nowParticle = particles[i]; //今回計算する粒子
        if (nowParticle.is_wall) continue;  // 速度の計算や位置の更新をしない粒子の場合スキップ
        let sum = createVector3(0, 0, 0); //足し合わせる変数
        for (let j = 0; j < particles.length; j++) { //他の粒子全てについて
            if(i == j){continue;} //自分自身だったらスキップ
            let nearParticle = particles[j];
            
            let diff = subVector3(nearParticle.position, nowParticle.position); //粒子距離
            const r2 = dotVector3(diff, diff); //粒子距離の２乗

            //粒子距離がhより小さい場合だけ計算する
            if ( r2 < h2 ) {
                const r = Math.sqrt(r2); //粒子距離
                const c = h - r;
                const n = ((nearParticle.pressure /*-*/+ nowParticle.pressure) / (2 * nearParticle.density)) * Math.pow(c,2) / r;
                sum = addVector3(sum, multiplyScalarVector3(diff, n));
            }
        }

        if (!terms[i]) terms[i] = {};
        terms[i].pressureTerm = multiplyScalarVector3(sum, (-1/*/nowParticle.pressure*/) * pressureCoef);  // 圧力項が求まった
    }
}

/**
 * 粒子の粘性項計算
 * @param {{position:{x:number, y:number, z:number}, velocity:{x:number, y:number, z:number}, force:{x:number, y:number, z:number}, density:number, pressure:number}[]} particles - 粒子のリスト
 */
function calcViscosityTerm(particles) {
    const h2 = h*h; //事前にhの二乗を計算しておく
    for (let i = 0; i < particles.length; i++) { //一つづつ粒子の密度を計算

        reportProgress(reportHeader + "粘性項を計算中..." + i + "/" + particles.length);

        let nowParticle = particles[i]; //今回計算する粒子
        if (nowParticle.is_wall) continue;  // 速度の計算や位置の更新をしない粒子の場合スキップ
        let sum = createVector3(0, 0, 0); //足し合わせる変数
        for (let j = 0; j < particles.length; j++) { //他の粒子全てについて
            if(i == j){continue;} //自分自身だったらスキップ
            let nearParticle = particles[j];
            
            let diff = subVector3(nearParticle.position, nowParticle.position); //粒子距離
            const r2 = dotVector3(diff, diff); //粒子距離の２乗

            //粒子距離がhより小さい場合だけ計算する
            if ( r2 < h2 ) {
                const r = Math.sqrt(r2); //粒子距離
                const c = h - r;
                const n = c / nearParticle.density;
                sum = addVector3(sum, multiplyScalarVector3(subVector3(nearParticle.velocity, nowParticle.velocity), n));
            }
        }

        if (!terms[i]) terms[i] = {};
        terms[i].viscosityTerm = multiplyScalarVector3(sum, viscosityCoef);  // 粘性項が求まった
    }
}


/**
 * 粒子のリスト
 * @type {{
 *      position:{x:number, y:number, z:number},
 *      velocity:{x:number, y:number, z:number},
 *      force:{x:number, y:number, z:number},
 *      acceleration:{x:number, y:number, z:number},
 *      density:number,
 *      pressure:number,
 *      is_wall:boolean
 * }[]}
 */
const _particles = [];


const deltaTime = 0.1;  // （秒）

let reportHeader = "";

function tick() {
    console.info("calcDensity");
    calcDensity(_particles);
    console.info("calcPressure");
    calcPressure(_particles);
    console.info("calcPressureTerm");
    calcPressureTerm(_particles);
    console.info("calcViscosityTerm");
    calcViscosityTerm(_particles);

    for (let i = 0; i < _particles.length; i++) {
        const nowParticle = _particles[i];
        if (nowParticle.is_wall) continue;
        const a = addVector3(addVector3(terms[i].pressureTerm, terms[i].viscosityTerm), g);
        const v = addVector3(nowParticle.velocity, multiplyScalarVector3(addVector3(nowParticle.acceleration, a), 0.5 * deltaTime));
        const deltaPosition = addVector3(multiplyScalarVector3(v, deltaTime), multiplyScalarVector3(a, 0.5 * deltaTime * deltaTime));
        nowParticle.acceleration = a;
        nowParticle.velocity = v;
        nowParticle.position = addVector3(nowParticle.position, deltaPosition);
    }
}

/**
 * シミュレーションを開始する関数
 * @param {number} simulateSeconds - シミュレーションを行う秒数
 */
function start(simulateSeconds) {
    // 初めに実行する処理
    makeWall(_particles);
    addParticles(_particles);

    self.postMessage({type:"result", content:_particles, time:0});

    for (let i = 0; i < simulateSeconds / deltaTime; i++) {
        const time = (i + 1) * deltaTime;
        reportHeader = time + "秒目/" + simulateSeconds + "秒 ";
        tick();
        self.postMessage({type:"result", content:_particles, time:time});
    }

    reportProgress("シミュレーション終了");
}
start(100);



/* 関数定義部分 */

/**
 * 粒子のオブジェクトを作成する関数
 * @param {{x:number, y:number, z:number}} position - 粒子の位置
 * @param {{x:number, y:number, z:number}} velocity - 粒子の速度
 * @param {{x:number, y:number, z:number}} force - 粒子にかかっている力
 * @param {{x:number, y:number, z:number}} acceleration - 粒子の加速度
 * @param {number} density - 粒子の密度
 * @param {number} pressure - 粒子の圧力
 * @param {boolean} is_wall - 粒子が壁かどうか
 * @returns {{position:{x:number, y:number, z:number}, velocity:{x:number, y:number, z:number}, force:{x:number, y:number, z:number}, acceleration:{x:number, y:number, z:number}, density:number, pressure:number, is_wall:boolean}} - 粒子のオブジェクト
 */
function createParticle(position={x:0,y:0,z:0}, is_wall=false, velocity={x:0,y:0,z:0}, force={x:0,y:0,z:0}, acceleration={x:0,y:0,z:0}, density=0, pressure=0) {
    return {
        position: {x: position.x, y: position.y, z: position.z},
        velocity: {x: velocity.x, y:velocity.y, z: velocity.z},
        force: {x: force.x, y: force.y, z: force.z},
        acceleration: {x: acceleration.x, y: acceleration.y, z: acceleration.z},
        density: density,
        pressure: pressure,
        is_wall: is_wall
    };
}

/**
 * ベクトルを表すオブジェクトを返す関数
 * @param {number} x - x成分
 * @param {number} y - y成分
 * @param {number} z - z成分
 * @returns - ベクトルを表すオブジェクト
 */
function createVector3(x = 0, y = 0, z = 0) {
    return {x: x, y: y, z: z};
}

/**
 * ベクトルについて a + b を計算する関数
 * @param {{x:number, y:number, z:number}} a - 足すベクトル
 * @param {{x:number, y:number, z:number}} b - 足されるベクトル
 * @returns {{x:number, y:number, z:number}} - ベクトルの和
 */
function addVector3(a, b) {
    return createVector3(a.x + b.x, a.y + b.y, a.z + b.z);
}

/**
 * ベクトルについて a - b を計算する関数
 * @param {{x:number, y:number, z:number}} a - 引かれるベクトル
 * @param {{x:number, y:number, z:number}} b - 引くベクトル
 * @returns {{x:number, y:number, z:number}} - ベクトルの差
 */
function subVector3(a, b) {
    return createVector3(a.x - b.x, a.y - b.y, a.z - b.z);
}

/**
 * ベクトルについて aのn倍 を計算する関数
 * @param {{x:number, y:number, z:number}} a - ベクトル
 * @param {number} n - 実数（スカラー）
 * @returns {{x:number, y:number, z:number}} - 実数倍したベクトル
 */
function multiplyScalarVector3(a, n) {
    return createVector3(n * a.x, n * a.y, n * a.z);
}

/**
 * ベクトルについて a ・ b （内積） を計算する関数
 * @param {{x:number, y:number, z:number}} a - ベクトル
 * @param {{x:number, y:number, z:number}} b - ベクトル
 * @returns {number} - ベクトルの内積
 */
function dotVector3(a, b) {
    return a.x * b.x + a.y * b.y + a.z * b.z;
}



/**
 * クライエントへ進捗を報告する関数
 * @param {string} text - 報告内容
 */
function reportProgress(text) {
    self.postMessage({type:"progress", content:text});
}