# Model providers

[Deep Agents Code](/oss/python/deepagents/code/overview)

# Model providers
Copy page

Configure any LangChain-compatible model provider for Deep Agents CodeCopy page
> 

## Documentation Index

Fetch the complete documentation index at: [https://docs.langchain.com/llms.txt](https://docs.langchain.com/llms.txt)

Use this file to discover all available pages before exploring further.Deep Agents Code supports any [chat model provider compatible with LangChain](/oss/python/integrations/chat), unlocking use for virtually any LLM that supports tool calling. Any service that exposes an OpenAI-compatible or Anthropic-compatible API also works out of the box—see [Compatible APIs](/oss/python/deepagents/code/configuration#compatible-apis).

## [​](#quickstart)Quickstart

Deep Agents Code integrates automatically with the [following model providers](#provider-reference): no extra configuration needed beyond installing the relevant provider package.

{' ' * (self.list_depth - 1)}- 
**Install provider packages**
Each model provider requires installing its corresponding LangChain integration package. These are available as optional extras when installing Deep Agents Code, done intentionally to keep the application lightweight:

```
# Quick install with chosen providers
# OpenAI, Anthropic, and Gemini are included by default
DEEPAGENTS_EXTRAS="baseten,groq" curl -LsSf https://langch.in/dcode | bash

# Or install directly with uv
uv tool install 'deepagents-code[baseten,groq]'

# Add additional packages at a later date
uv tool install deepagents-code --with langchain-ollama

# All providers
uv tool install 'deepagents-code[anthropic,baseten,bedrock,cohere,deepseek,fireworks,google-genai,groq,huggingface,ibm,litellm,mistralai,nvidia,ollama,openai,openrouter,perplexity,vertexai,xai]'

```

{' ' * (self.list_depth - 1)}- 
**Set credentials**
Store API keys in `~/.deepagents/.env` so they’re available across all projects, or export them in your shell:

{' ' * (self.list_depth - 1)}- OpenAI
{' ' * (self.list_depth - 1)}- Anthropic
{' ' * (self.list_depth - 1)}- Google
{' ' * (self.list_depth - 1)}- OtherAdd permanentlyAdd for current session
```
mkdir -p ~/.deepagents
echo 'OPENAI_API_KEY=your-api-key' >> ~/.deepagents/.env

```
Add permanentlyAdd for current session
```
mkdir -p ~/.deepagents
echo 'ANTHROPIC_API_KEY=your-api-key' >> ~/.deepagents/.env

```
Add permanentlyAdd for current session
```
mkdir -p ~/.deepagents
echo 'GOOGLE_API_KEY=your-api-key' >> ~/.deepagents/.env

```
Deep Agents Code works with any LLM that supports tool calling. See [Provider reference](#provider-reference) for the full list of supported providers and their required environment variables.
To configure model parameters, see [Model parameters](/oss/python/deepagents/code/providers#model-parameters).
You can also scope credentials to Deep Agents Code with the [`DEEPAGENTS_CODE_` prefix](/oss/python/deepagents/code/configuration#deepagents_code_-prefix).

## [​](#provider-reference)Provider reference

Using a provider not listed here? See [Arbitrary providers](/oss/python/deepagents/code/configuration#arbitrary-providers): any LangChain-compatible provider can be used in Deep Agents Code with additional setup.

| 
 | Provider | Package | Credential env var | Model profiles
 | OpenAI | [`langchain-openai`](/oss/python/integrations/chat/openai) | `OPENAI_API_KEY` | ✅
 | Azure OpenAI | [`langchain-openai`](/oss/python/integrations/chat/azure_chat_openai) | `AZURE_OPENAI_API_KEY` | ✅
 | Anthropic | [`langchain-anthropic`](/oss/python/integrations/chat/anthropic) | `ANTHROPIC_API_KEY` | ✅
 | Google Gemini API | [`langchain-google-genai`](/oss/python/integrations/chat/google_generative_ai) | `GOOGLE_API_KEY` | ✅
 | Google Vertex AI | [`langchain-google-genai`](/oss/python/integrations/chat/google_generative_ai#credentials) | `GOOGLE_CLOUD_PROJECT` | ✅
 | Baseten | [`langchain-baseten`](https://github.com/basetenlabs/langchain-baseten) | `BASETEN_API_KEY` | ✅
 | AWS Bedrock | [`langchain-aws`](/oss/python/integrations/chat/bedrock) | `AWS_ACCESS_KEY_ID`, `AWS_SECRET_ACCESS_KEY` | ✅
 | AWS Bedrock Converse | [`langchain-aws`](/oss/python/integrations/chat/bedrock) | `AWS_ACCESS_KEY_ID`, `AWS_SECRET_ACCESS_KEY` | ✅
 | Hugging Face | [`langchain-huggingface`](/oss/python/integrations/chat/huggingface) | `HUGGINGFACEHUB_API_TOKEN` | ✅
 | Ollama | [`langchain-ollama`](/oss/python/integrations/chat/ollama) | `OLLAMA_API_KEY` (cloud only; optional) | ❌
 | Groq | [`langchain-groq`](/oss/python/integrations/chat/groq) | `GROQ_API_KEY` | ✅
 | Cohere | [`langchain-cohere`](/oss/python/integrations/chat/cohere) | `COHERE_API_KEY` | ❌
 | Fireworks | [`langchain-fireworks`](/oss/python/integrations/chat/fireworks) | `FIREWORKS_API_KEY` | ✅
 | Together | [`langchain-together`](/oss/python/integrations/chat/together) | `TOGETHER_API_KEY` | ❌
 | Mistral AI | [`langchain-mistralai`](/oss/python/integrations/chat/mistralai) | `MISTRAL_API_KEY` | ✅
 | DeepSeek | [`langchain-deepseek`](/oss/python/integrations/chat/deepseek) | `DEEPSEEK_API_KEY` | ✅
 | IBM (watsonx.ai) | [`langchain-ibm`](/oss/python/integrations/chat/ibm_watsonx) | `WATSONX_APIKEY` | ❌
 | Nvidia | [`langchain-nvidia-ai-endpoints`](/oss/python/integrations/chat/nvidia_ai_endpoints) | `NVIDIA_API_KEY` | ✅
 | xAI | [`langchain-xai`](/oss/python/integrations/chat/xai) | `XAI_API_KEY` | ✅
 | Perplexity | [`langchain-perplexity`](/oss/python/integrations/chat/perplexity) | `PERPLEXITY_API_KEY` (or `PPLX_API_KEY`) | ✅
 | OpenRouter | [`langchain-openrouter`](/oss/python/integrations/chat/openrouter) | `OPENROUTER_API_KEY` | ✅
 | LiteLLM | [`langchain-litellm`](/oss/python/integrations/chat/litellm) | Per-provider (see [docs](https://docs.litellm.ai/)) | ❌
You can scope any credential to Deep Agents Code by adding a `DEEPAGENTS_CODE_` prefix. For example, `DEEPAGENTS_CODE_OPENAI_API_KEY` takes priority over `OPENAI_API_KEY` within Deep Agents Code without affecting other tools. See [`DEEPAGENTS_CODE_` prefix](/oss/python/deepagents/code/configuration#deepagents_code_-prefix) for details.
A **[model profile](/oss/python/langchain/models#model-profiles)** is a bundle of metadata (model name, default parameters, capabilities, etc.) that ships with a provider package, largely powered by the [models.dev](https://models.dev/) project.Providers that include model profiles have their models listed automatically in the interactive `/model` switcher, subject to the [filtering criteria](#which-models-appear-in-the-switcher) (notably, `tool_calling` must be enabled). Providers without model profiles require you to specify the model name directly or add models via `config.toml`.

### [​](#model-routers-and-proxies)Model routers and proxies

Model routers like [OpenRouter](https://openrouter.ai/) and [LiteLLM](https://docs.litellm.ai/) provide access to models from multiple providers through a single endpoint.
Use the dedicated integration packages for these services:

| 
 | Router | Package | Config
 | OpenRouter | [`langchain-openrouter`](/oss/python/integrations/chat/openrouter) | `openrouter:<model>` (built-in, see [Provider reference](#provider-reference))
 | LiteLLM | [`langchain-litellm`](/oss/python/integrations/chat/litellm) | `litellm:<model>` (built-in, see [Provider reference](#provider-reference))
**OpenRouter** is a built-in provider—install the package and use it directly:

```
uv tool install 'deepagents-code[openrouter]'

```

**LiteLLM** is also a built-in provider:

```
uv tool install 'deepagents-code[litellm]'

```

## [​](#switch-models)Switch models

To switch models in Deep Agents Code, either:

{' ' * (self.list_depth - 1)}- 
**Use the interactive model switcher** with the `/model` command. This displays available models sourced from each installed LangChain provider package’s [model profiles](/oss/python/langchain/models#model-profiles).
Not all models appear here. If yours is missing, pass the model name directly (e.g. `/model gpt-5.5`). See [Which models appear in the switcher](#which-models-appear-in-the-switcher) for details.

{' ' * (self.list_depth - 1)}- 
**Specify a model name directly** as an argument, e.g. `/model gpt-5.5`. You can use any model supported by the chosen provider, regardless of whether it appears in the list from option 1. The model name will be passed to the API request.

{' ' * (self.list_depth - 1)}- 
**Specify the model at launch** via `--model`, e.g.

```
dcode --model openai:gpt-5.5

```

Model resolution orderWhen Deep Agents Code launches, it resolves which model to use in the following order:

{' ' * (self.list_depth - 1)}- **`--model` flag** always wins when provided.

{' ' * (self.list_depth - 1)}- **`[models].default`** in `~/.deepagents/config.toml`—the user’s intentional long-term preference.

{' ' * (self.list_depth - 1)}- **`[models].recent`** in `~/.deepagents/config.toml`—the last model switched to via `/model`. Written automatically; never overwrites `[models].default`.

{' ' * (self.list_depth - 1)}- **Environment auto-detection**: falls back to the first available startup credential, checked in order: `OPENAI_API_KEY`, `ANTHROPIC_API_KEY`, `GOOGLE_API_KEY`, `GOOGLE_CLOUD_PROJECT` (Vertex AI).
This startup fallback intentionally checks only those four credentials. Other supported providers (for example, Groq) are still available via `--model`, `/model`, and saved defaults (`[models].default` / `[models].recent`).

### [​](#which-models-appear-in-the-switcher)Which models appear in the switcher

The `/model` selector dynamically builds its list from installed provider packages. Expand below for the full criteria and troubleshooting.

How the switcher builds its model listThe interactive `/model` selector builds its list dynamically—it is not a hardcoded list baked into Deep Agents Code. A model appears in the switcher when **all** of the following are true:

{' ' * (self.list_depth - 1)}- 
**The provider package is installed.** Each provider (e.g. `langchain-anthropic`, `langchain-openai`) must be installed alongside `deepagents-code`—either as an [install extra](/oss/python/deepagents/code/providers#quickstart) (e.g. `uv tool install 'deepagents-code[ollama]'`) or added later with `uv tool install deepagents-code --with <package>`. If a package is missing, its entire provider section is absent from the switcher.

{' ' * (self.list_depth - 1)}- 
**The model has a profile with `tool_calling` enabled.** Deep Agents Code requires tool-calling support, so models without `tool_calling: true` in their profile are excluded. This is the most common reason a model is missing from the list. For providers that don’t bundle profiles (see the [Provider reference](#provider-reference) table), you can define one in `config.toml`:

```
[models.providers.ollama.profile."qwen3:4b"]
tool_calling = true
max_input_tokens = 32768
max_output_tokens = 8192

```

This is not strictly required for the model to appear in the switcher — adding it to the [`models` list](/oss/python/deepagents/code/configuration#adding-models-to-the-interactive-switcher) also works and is simpler. A profile is useful when you want Deep Agents Code to know the model’s context window and capabilities for features like auto-summarization. See [Profile overrides](/oss/python/deepagents/code/configuration#profile-overrides-advanced) for all overridable fields.

{' ' * (self.list_depth - 1)}- 
**The model accepts and produces text.** Models whose profile explicitly sets `text_inputs` or `text_outputs` to `false` (e.g. embedding or image-generation models) are excluded.

Models defined in `config.toml` under [`[models.providers.<name>].models`](/oss/python/deepagents/code/configuration#adding-models-to-the-interactive-switcher) bypass the profile filter—they always appear in the switcher regardless of profile metadata. This is the recommended way to add models that are missing from the list.Credential status does **not** affect whether a model is listed. The switcher shows all qualifying models and displays a credential indicator next to each provider header: a checkmark for confirmed credentials, a warning for missing credentials, or a question mark when credential status is unknown. You can still select a model with missing credentials—the provider will report an authentication error at request time.

#### [​](#troubleshooting-missing-models)Troubleshooting missing models

| 
 | Symptom | Likely cause | Fix
 | Entire provider missing from switcher | Provider package not installed | Install the package (e.g. `uv tool install deepagents-code --with langchain-groq`)
 | Provider shown but specific model missing | Model profile has `tool_calling: false` or no profile exists | Add the model to `[models.providers.<name>].models` in `config.toml`, or use `/model <provider>:<model>` directly
 | Provider shows ⚠ “missing credentials” | API key env var not set | Set the credential env var from the [Provider reference](#provider-reference) table
 | Provider shows ? “credentials unknown” | Provider uses non-standard auth that Deep Agents Code can’t verify | Credentials may still work—try switching to the model. If auth fails, check the provider’s docs

### [​](#set-a-default-model)Set a default model

You can set a persistent default model that will be used for all future CLI launches:

{' ' * (self.list_depth - 1)}- 
**Via model selector:** Open `/model`, navigate to the desired model, and press `Ctrl+S` to pin it as the default. Pressing `Ctrl+S` again on the current default clears it.

{' ' * (self.list_depth - 1)}- 
**Via command:** `/model --default provider:model` (e.g., `/model --default anthropic:claude-opus-4-7`)

{' ' * (self.list_depth - 1)}- 
**Via config file:** Set `[models].default` in `~/.deepagents/config.toml` (see [Configuration](/oss/python/deepagents/code/configuration)).

{' ' * (self.list_depth - 1)}- 
**From the shell:**

```
dcode --default-model anthropic:claude-opus-4-7

```

To view the current default:

```
dcode --default-model

```

To clear the default:

{' ' * (self.list_depth - 1)}- 
**From the shell:**

```
dcode --clear-default-model

```

{' ' * (self.list_depth - 1)}- 
**Via command:** `/model --default --clear`

{' ' * (self.list_depth - 1)}- 
**Via model selector:** Press `Ctrl+S` on the currently pinned default model.

Without a default, Deep Agents Code will default to the most recently used model.

### [​](#model-parameters)Model parameters

Pass extra constructor kwargs to the model—sampling controls, reasoning/thinking budgets, context window sizes, request timeouts, and anything else the underlying chat-model class accepts. Three places to set them, in priority order (highest first):

{' ' * (self.list_depth - 1)}- 
**One-off at launch with `--model-params`.** JSON string, session-only:

```
# OpenAI reasoning effort
dcode --model openai:gpt-5.5 --model-params '{"reasoning": {"effort": "high"}}'

# Anthropic extended thinking
dcode --model anthropic:claude-opus-4-7 --model-params '{"thinking": {"type": "enabled", "budget_tokens": 10000}, "max_tokens": 16000}'

```

{' ' * (self.list_depth - 1)}- 
**Mid-session via `/model --model-params`.** Same JSON syntax—swaps params (and optionally the model) without restarting:

```
/model --model-params '{"temperature": 0.7}' anthropic:claude-opus-4-7
/model --model-params '{"num_ctx": 16384}' # opens selector, applies params to choice

```

{' ' * (self.list_depth - 1)}- 
**Persistent in `config.toml`.** Provider-level defaults (with optional per-model sub-tables) that apply on every launch:

```
[models.providers.anthropic.params]
thinking = { type = "enabled", budget_tokens = 10000 }
max_tokens = 16000

[models.providers.openai.params]
reasoning = { effort = "high", summary = "auto" }
output_version = "responses/v1"

[models.providers.ollama.params]
num_ctx = 16384
temperature = 0

# Per-model override—wins over provider-level keys
[models.providers.ollama.params."qwen3:4b"]
temperature = 0.5

```

CLI flags override config-file `params` and are session-only (mid-session changes are not persisted). Per-model sub-tables in `config.toml` override provider-level keys (shallow merge—see [Model constructor params](/oss/python/deepagents/code/configuration#model-constructor-params) for full semantics). `--model-params` cannot be combined with `--default`.
Any kwarg accepted by the underlying chat-model constructor is valid. Refer to the provider’s reference docs for the full list—e.g. [`ChatAnthropic`](https://reference.langchain.com/python/langchain-anthropic/langchain_anthropic/chat_models/ChatAnthropic), [`ChatOpenAI`](https://reference.langchain.com/python/langchain-openai/langchain_openai/chat_models/base/ChatOpenAI), [`ChatOllama`](https://reference.langchain.com/python/langchain-ollama/langchain_ollama/chat_models/ChatOllama). Unknown kwargs are forwarded to the upstream API request, so newly released parameters work without a CLI update.
Don’t put credentials (`api_key`) in `params`—use [`api_key_env`](/oss/python/deepagents/code/configuration#provider-configuration) to point at an environment variable instead.
To override fields on the model’s runtime *profile* (`max_input_tokens`, `tool_calling`, capability flags)—distinct from constructor params—see [Profile overrides](/oss/python/deepagents/code/configuration#profile-overrides-advanced).

## [​](#advanced-configuration)Advanced configuration

For detailed configuration of provider params, profile overrides, custom base URLs, compatible APIs, arbitrary providers, and lifecycle hooks, see [Configuration](/oss/python/deepagents/code/configuration).

---

[Connect these docs](/use-these-docs) to Claude, VSCode, and more via MCP for real-time answers.[Edit this page on GitHub](https://github.com/langchain-ai/docs/edit/main/src/oss/deepagents/code/providers.mdx) or [file an issue](https://github.com/langchain-ai/docs/issues/new/choose).

Was this page helpful?YesNo[Use subagents in Deep Agents CodePrevious](/oss/python/deepagents/code/subagents)[ConfigurationNext](/oss/python/deepagents/code/configuration)⌘I[Docs by LangChain home page
![light logo](https://mintcdn.com/langchain-5e9cc07a/nQm-sjd_MByLhgeW/images/brand/langchain-docs-dark-blue.png?fit=max&auto=format&n=nQm-sjd_MByLhgeW&q=85&s=5babf1a1962208fd7eed942fa2432ecb)
![dark logo](https://mintcdn.com/langchain-5e9cc07a/nQm-sjd_MByLhgeW/images/brand/langchain-docs-light-blue.png?fit=max&auto=format&n=nQm-sjd_MByLhgeW&q=85&s=0bcd2a1f2599ed228bcedf0f535b45b1)](/)[github](https://github.com/langchain-ai)[x](https://x.com/LangChain)[linkedin](https://www.linkedin.com/company/langchain)[youtube](https://www.youtube.com/@LangChain)

Resources[Forum](https://forum.langchain.com/)[Changelog](https://changelog.langchain.com/)[LangChain Academy](https://academy.langchain.com/)[Contact Sales](https://www.langchain.com/contact-sales)

Company[Home](https://langchain.com/)[Trust Center](https://trust.langchain.com/)[Careers](https://langchain.com/careers)[Blog](https://blog.langchain.com/)[github](https://github.com/langchain-ai)[x](https://x.com/LangChain)[linkedin](https://www.linkedin.com/company/langchain)[youtube](https://www.youtube.com/@LangChain)
